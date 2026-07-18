import pandas as pd
from openai import OpenAI
import pronto
from celltyper_data import data_file
from langgraph.graph import START, StateGraph, END
from cell_marker_graph import construct_marker_graph
import hashlib
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from dotenv import load_dotenv
from typing import Optional
from igraph import Graph
from tqdm import tqdm
from geneRetriever import scene_gene_desc_summary, scene_gene_function, scene_cell_relatedTo_gene, scene_pathway, scene_protein_location, p_p_interaction_info, unknown_genes_info, route_graph_construction, route_unknown_genes, normalize_arabidopsis_markers_str, get_gene_db
from typing import TypedDict


LLM_SELF_CONSISTENCY_SAMPLES = 5


load_dotenv()







class CellTyperState(TypedDict):
  markers: Optional[str] = None
  species: Optional[str] = None
  tissue: Optional[str] = None
  condition: Optional[str] = None
  db1: Optional[pd.DataFrame] = None
  db2: Optional[pd.DataFrame] = None
  db3: Optional[pd.DataFrame] = None
  db5: Optional[pd.DataFrame] = None
  db6: Optional[pd.DataFrame] = None
  db8: Optional[pd.DataFrame] = None
  geneDB: Optional[pd.DataFrame] = None
  ont: Optional[pronto.Ontology] = None
  marker_cell_graph: Optional[Graph] = None
  celltype_candidates: Optional[list] = None   
  fin_celltype: Optional[str] = None
  pbar: Optional[tqdm] = None
  gene_desc_summary_context: Optional[str] = None
  genes_ids: Optional[dict] = None
  gene_func_summary_context: Optional[str] = None
  gene_cell_summary_context: Optional[str] = None
  gene_pathway_context: Optional[str] = None
  gene_protein_loc_context: Optional[str] = None
  ppi_context: Optional[str]=None
  genetriever_flag: Optional[bool] = None
  unknown_genes: Optional[list] = None
  llm_model: Optional[str] = None
  graph_provided: Optional[bool] = None
  llm_sample_texts: Optional[list] = None
  llm_parsed_celltypes: Optional[list] = None
  llm_consistency_counts: Optional[dict] = None
  llm_consistency_scores: Optional[dict] = None
  fin_celltype_display: Optional[str] = None
  cluster_id: Optional[str] = None
  progress_disable: Optional[bool] = None

def prepare(state):
  cluster_id = state.get("cluster_id")
  if cluster_id is not None and str(cluster_id).strip():
    desc = f"Cluster {cluster_id}"
  else:
    desc = "Progress"
  disable = bool(state.get("progress_disable"))
  pbar = tqdm(
    total=100,
    desc=desc,
    bar_format="{l_bar}{bar}| {percentage:3.0f}%",
    disable=disable,
    leave=not disable,
  )

  pbar.update(5)

  markers = state.get("markers")
  species = state.get("species")
  geneDB = state.get("geneDB")
  if geneDB is None and species:
    geneDB = get_gene_db(species)
  if markers and species and str(species).lower() == "arabidopsis":
    markers = normalize_arabidopsis_markers_str(markers, db8=state.get("db8"))
  
  return {**state, "pbar":pbar, "markers":markers, "geneDB":geneDB}

def find_celltypes_in_tissue(state):
  species = state["species"]
  tissue = state["tissue"]
  celltypes = []
  
  if species == "human":
    species = "homo_sapiens"
    with open(data_file("human_celltype.json"), "r", encoding="utf-8") as f:
      data = json.load(f)
    celltypes = data[tissue]
    
  if species == "mouse":
    species = "mus_musculus"
    with open(data_file("musmusculus_celltype.json"), "r", encoding="utf-8") as f:
      data = json.load(f)
    celltypes = data[tissue]
  
  if species == "arabidopsis":
    with open(data_file("arabidopsis_celltype.json"), "r", encoding="utf-8") as f:
      data = json.load(f)
    
    if tissue == "shoot":
      celltypes = data["leaf"] | data["fruit"] | data["flower"] | data["shoot apex"] | data["stem"]
    elif tissue == "seedling":
      celltypes = data["leaf"] | data["root"] | data["shoot apex"] | data["stem"]
    else:
      celltypes = data[tissue]
   
  pbar = state["pbar"]
  pbar.update(5)
  return {**state, "celltype_candidates":celltypes, "pbar":pbar}
  


def _parse_cell_types_from_llm_text(text):
  if not text or not str(text).strip():
    return []
  for marker in ("Cell types:", "**Cell types**:"):
    pos = text.find(marker)
    if pos < 0:
      continue
    tail = text[pos + len(marker) :].split("\n", 1)[0]
    parts = [item.strip("** `").strip() for item in tail.split("&")]
    out = [p for p in parts if p]
    if out:
      return out
  m = re.search(
    r"(?:\*\*)?Cell types?(?:\*\*)?\s*:\s*([^\n]+)",
    text,
    flags=re.IGNORECASE,
  )
  if m:
    tail = m.group(1)
    parts = [item.strip("** `").strip() for item in tail.split("&")]
    return [p for p in parts if p]
  return []


def _self_consistency_prompt_cache_key(prompt: str) -> str:
  return hashlib.sha256(prompt.encode("utf-8")).hexdigest()


def _is_reasoning_model(model: str) -> bool:
  m = (model or "").strip().lower()
  if not m:
    return False
  return (
    m.startswith("gpt-5")
    or m.startswith("o1")
    or m.startswith("o3")
    or m.startswith("o4")
  )


def _openai_text_from_messages(
  llm,
  model: str,
  messages,
  *,
  temperature: float = 0.8,
  prompt_cache_key: Optional[str] = None,
) -> str:
  model = (model or "").strip()
  if not model:
    return ""

  if _is_reasoning_model(model):
    kwargs = {
      "model": model,
      "input": messages,
      "reasoning": {"effort": "low"},
    }
    if prompt_cache_key:
      kwargs["prompt_cache_key"] = prompt_cache_key
    try:
      response = llm.responses.create(**kwargs)
      text = getattr(response, "output_text", None) or ""
      if text.strip():
        return text
    except Exception:
      pass

  chat_kwargs = {
    "model": model,
    "messages": messages,
    "temperature": temperature,
  }
  if prompt_cache_key:
    chat_kwargs["extra_body"] = {"prompt_cache_key": prompt_cache_key}
  try:
    completion = llm.chat.completions.create(**chat_kwargs)
    text = completion.choices[0].message.content or ""
    if text.strip() or not _is_reasoning_model(model):
      return text
  except Exception:
    if not _is_reasoning_model(model):
      return ""

  if _is_reasoning_model(model):
    try:
      response = llm.responses.create(
        model=model,
        input=messages,
        reasoning={"effort": "low"},
      )
      return getattr(response, "output_text", None) or ""
    except Exception:
      return ""
  return ""


def _run_self_consistency_sample(llm, llm_model, messages, prompt_cache_key):
  return _openai_text_from_messages(
    llm,
    llm_model,
    messages,
    temperature=0.8,
    prompt_cache_key=prompt_cache_key,
  )


def _collect_self_consistency_samples(llm, llm_model, messages, prompt_cache_key):
  n = LLM_SELF_CONSISTENCY_SAMPLES
  sample_texts = [""] * n
  if n <= 0:
    return sample_texts

  sample_texts[0] = _run_self_consistency_sample(
    llm, llm_model, messages, prompt_cache_key
  )
  if n == 1:
    return sample_texts

  with ThreadPoolExecutor(max_workers=n - 1) as executor:
    futures = {
      executor.submit(
        _run_self_consistency_sample, llm, llm_model, messages, prompt_cache_key
      ): idx
      for idx in range(1, n)
    }
    for future in as_completed(futures):
      sample_texts[futures[future]] = future.result()

  return sample_texts


def select_candidates_from_llm(state):
  llm = OpenAI()
  cell_candidates = state["celltype_candidates"]
  tissue = state["tissue"]
  species = state["species"]
  markers = state["markers"]
  condition = state["condition"]
  raw_model = state.get("llm_model")
  if raw_model is None:
    llm_model = ""
  else:
    llm_model = str(raw_model).strip()
  if not llm_model:
    llm_model = (
      os.environ.get("CELLTYPER_LLM_MODEL")
      or os.environ.get("OPENAI_MODEL")
      or "gpt-5.2"
    ).strip()
  gene_desc_summary_context = state["gene_desc_summary_context"]
  gene_func_summary_context = state["gene_func_summary_context"]
  gene_cell_summary_context = state["gene_cell_summary_context"]
  gene_pathway_context = state["gene_pathway_context"]
  gene_protein_loc_context = state["gene_protein_loc_context"]
  ppi_context = state.get("ppi_context") or ""
  
  
  prompt = f"""You are an intelligent cell annotation agent that identifies the most likely cell type candidate in the given marker genes and conditions. You should think step by step through the given marker genes, their function, and their combinations.

  
  You will be given total cell type candidates and you should choose the most likely cell type considering the given marker genes.

  
  You can use the given following context:
  1. Gene Description & Summary
  2. Gene Function
  3. Some Cell Types related to Gene
  4. Pathway related to Gene
  5. Subcellular Location of Gene's Product
  6. Protein-Protein Interaction
  
  You should not pick the cell type with the reason that lots of certain cell type can be observed in the given tissue.
  
  The last sentence from you should be Cell types: [Cell name] (ex. Cell types: Acinar cell)
  
  You should strictly distinguish the given cell type candidates with given markers
  
  The final cell type you pick should be supported by the stronger signal of marker genes than other cell type supported by the weaker signal of marker genes.
  
  
  
  [Cell type candidates from {tissue}]
  {json.dumps(cell_candidates, indent=2)}
  
  [Marker genes]
  {markers}
  
  [Condition]
  species: {species}
  Additional: {condition}
  
  """
  
  
  if species == "arabidopsis":
    prompt += "If the tissue is shoot, you should consider not only leaf but also stem, shoot apex, fruit, flower and so on.\n"
    prompt += "**Note that you should consider the cell type information related to the gene more than the biological process or function of gene.**\n"
    prompt += "**Note that you are given cell candidates with a tree structure, so you should start from the top-level categories and reason in a top-down elimination manner. You should stop at the most appropriate level of cell type specificity supported by the evidence.**"
    
  prompt += gene_desc_summary_context + "\n" + gene_func_summary_context  + "\n" + gene_cell_summary_context + "\n" + gene_pathway_context + gene_protein_loc_context + "\n" + ppi_context + "\n"
  
  
  messages = [{"role": "user", "content": prompt}]
  prompt_cache_key = _self_consistency_prompt_cache_key(prompt)
  sample_texts = _collect_self_consistency_samples(
    llm, llm_model, messages, prompt_cache_key
  )

  rank = {}
  parsed_celltypes = []
  for text in sample_texts:
    labels = _parse_cell_types_from_llm_text(text)
    primary = labels[0] if labels else ""
    parsed_celltypes.append(primary)

    if primary:
      key = primary.lower()
      rank[key] = rank.get(key, 0) + 1

  
  label_map = {}
  for p in parsed_celltypes:
    if p and p.lower() not in label_map:
      label_map[p.lower()] = p
  for k in rank:
    if k not in label_map:
      label_map[k] = k

  if rank:
    max_votes = max(rank.values())
    winner_keys = sorted(
      [k for k, v in rank.items() if v == max_votes],
      key=lambda k: label_map.get(k, k).lower(),
    )
  else:
    winner_keys = []

  winner_labels = [label_map[k] for k in winner_keys]
  fin_celltype_display = " & ".join(winner_labels)
  fin_celltype = " & ".join(winner_keys)

  consistency_counts = {label_map[k]: rank[k] for k in rank}
  n = float(LLM_SELF_CONSISTENCY_SAMPLES)
  consistency_scores = {label: consistency_counts[label] / n for label in consistency_counts}

  
  pbar = state["pbar"]
  pbar.update(pbar.total - pbar.n)
  
  return {
    **state,
    "fin_celltype": fin_celltype,
    "fin_celltype_display": fin_celltype_display,
    "llm_sample_texts": sample_texts,
    "llm_parsed_celltypes": parsed_celltypes,
    "llm_consistency_counts": consistency_counts,
    "llm_consistency_scores": consistency_scores,
    "pbar": pbar,
  }
  


    




CellTyper_workflow = StateGraph(CellTyperState)

CellTyper_workflow.add_node("prepare", prepare)
CellTyper_workflow.add_node("Construct_Cell_Marker_Graph", construct_marker_graph)
CellTyper_workflow.add_node("find_celltypes_in_tissue", find_celltypes_in_tissue)
CellTyper_workflow.add_node("select_candidates_from_llm", select_candidates_from_llm)

CellTyper_workflow.add_node("scene_gene_desc_summary", scene_gene_desc_summary)
CellTyper_workflow.add_node("scene_gene_function", scene_gene_function)
CellTyper_workflow.add_node("scene_cell_relatedTo_gene", scene_cell_relatedTo_gene)
CellTyper_workflow.add_node("scene_pathway", scene_pathway)
CellTyper_workflow.add_node("scene_protein_location", scene_protein_location)
CellTyper_workflow.add_node("p_p_interaction_info", p_p_interaction_info)
CellTyper_workflow.add_node("unknown_genes_info", unknown_genes_info)

CellTyper_workflow.add_edge(START, "prepare")
CellTyper_workflow.add_edge("prepare", "find_celltypes_in_tissue")
CellTyper_workflow.add_edge("find_celltypes_in_tissue", "scene_gene_desc_summary")
CellTyper_workflow.add_edge("scene_gene_desc_summary", "scene_gene_function")
CellTyper_workflow.add_conditional_edges(
  "scene_gene_function",
  route_graph_construction,
  {
    "Construct_Cell_Marker_Graph": "Construct_Cell_Marker_Graph",
    "scene_cell_relatedTo_gene": "scene_cell_relatedTo_gene"
  }
)

CellTyper_workflow.add_edge("Construct_Cell_Marker_Graph", "scene_cell_relatedTo_gene")
CellTyper_workflow.add_edge("scene_cell_relatedTo_gene", "scene_pathway")
CellTyper_workflow.add_edge("scene_pathway", "scene_protein_location")
CellTyper_workflow.add_edge("scene_protein_location", "p_p_interaction_info")
CellTyper_workflow.add_conditional_edges(
  "p_p_interaction_info",
  route_unknown_genes,
  {
    "select_candidates_from_llm": "select_candidates_from_llm",
    "unknown_genes_info": "unknown_genes_info"
  }
)

CellTyper_workflow.add_edge("unknown_genes_info", "select_candidates_from_llm")

CellTyper_workflow.add_edge("select_candidates_from_llm", END)

app = CellTyper_workflow.compile()






  
  
  