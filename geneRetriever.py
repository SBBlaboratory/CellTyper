import re
import requests
import pandas as pd
import mygene
import io, contextlib
from typing import Dict, List, Optional, Set

from celltyper_data import data_file
from dotenv import load_dotenv
from openai import OpenAI
from xml.etree import ElementTree

load_dotenv()


_TAIR_ID_RE = re.compile(r"^AT([1-5]G|CG|MG)\d{5}(\.\d+)?$", re.IGNORECASE)
_TAIR_ID_SEARCH_RE = re.compile(r"AT(?:[1-5]G|CG|MG)\d{5}", re.IGNORECASE)

_ARABIDOPSIS_GENEDB_GENES: Optional[Set[str]] = None
_ARABIDOPSIS_DB8_SYMBOL_MAP: Optional[Dict[str, str]] = None
_GENE_DB_CACHE: Dict[str, pd.DataFrame] = {}
_NCBI_TAIR_CACHE: Dict[str, List[str]] = {}
_MAPPING_RESULT_CACHE: Dict[str, str] = {}


def _normalize_gene_db_species(species: str) -> str:
  key = str(species).lower().strip()
  if key == "arabidopsis thaliana":
    return "arabidopsis"
  return key


def get_gene_db(species: str) -> pd.DataFrame:
  key = _normalize_gene_db_species(species)
  if key in _GENE_DB_CACHE:
    return _GENE_DB_CACHE[key]

  paths = {
    "human": "human_geneDB.csv",
    "mouse": "mouse_geneDB.csv",
    "arabidopsis": "arabidopsis_geneDB.csv",
  }
  if key not in paths:
    raise ValueError(f"Unsupported species for geneDB: {species}")

  _GENE_DB_CACHE[key] = pd.read_csv(data_file(paths[key]))
  return _GENE_DB_CACHE[key]


def get_gene_db_from_state(state) -> pd.DataFrame:
  gene_db = state.get("geneDB")
  if gene_db is not None:
    return gene_db
  return get_gene_db(state["species"])


def is_arabidopsis_tair_id(gene) -> bool:
  return bool(_TAIR_ID_RE.match(str(gene).strip()))


def _arabidopsis_tair_base(gene: str) -> str:
  return str(gene).strip().upper().split(".")[0]


def _arabidopsis_gene_db_genes() -> Set[str]:
  global _ARABIDOPSIS_GENEDB_GENES
  if _ARABIDOPSIS_GENEDB_GENES is None:
    db = get_gene_db("arabidopsis")
    _ARABIDOPSIS_GENEDB_GENES = set(db["Gene"].astype(str).str.upper())
  return _ARABIDOPSIS_GENEDB_GENES


def _arabidopsis_db8_symbol_map(db8) -> Dict[str, str]:
  global _ARABIDOPSIS_DB8_SYMBOL_MAP
  if _ARABIDOPSIS_DB8_SYMBOL_MAP is not None:
    return _ARABIDOPSIS_DB8_SYMBOL_MAP
  if db8 is None or "name" not in db8.columns or "gene" not in db8.columns:
    return {}
  mp: Dict[str, str] = {}
  subset = db8.dropna(subset=["name", "gene"]).drop_duplicates(subset=["name"], keep="first")
  for name, gene in zip(subset["name"].astype(str), subset["gene"].astype(str)):
    key = name.strip().lower()
    if not key:
      continue
    gene_s = gene.strip()
    mp[key] = _arabidopsis_tair_base(gene_s) if is_arabidopsis_tair_id(gene_s) else gene_s
  _ARABIDOPSIS_DB8_SYMBOL_MAP = mp
  return mp


def _arabidopsis_tair_ids_from_ncbi(symbol: str) -> List[str]:
  key = symbol.strip().lower()
  if key in _NCBI_TAIR_CACHE:
    return _NCBI_TAIR_CACHE[key]

  esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
  params = {
    "db": "gene",
    "term": f"{symbol}[gene] AND Arabidopsis thaliana[orgn]",
    "retmode": "json",
  }
  uniq: List[str] = []
  try:
    res = requests.get(esearch_url, params=params, timeout=20).json()
    ids = res.get("esearchresult", {}).get("idlist", [])
    if ids:
      efetch = requests.get(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
        params={"db": "gene", "id": ids[0], "retmode": "xml"},
        timeout=20,
      )
      seen = set()
      for match in _TAIR_ID_SEARCH_RE.findall(efetch.text):
        tair = match.upper()
        if tair not in seen:
          seen.add(tair)
          uniq.append(tair)
  except Exception:
    pass

  _NCBI_TAIR_CACHE[key] = uniq
  return uniq


def _pick_preferred_tair(candidates: List[str], gene_db_genes: Set[str]) -> Optional[str]:
  for tair in candidates:
    if tair in gene_db_genes:
      return tair
  return candidates[0] if candidates else None


def map_arabidopsis_gene_to_tair(
  gene,
  *,
  db8=None,
  db8_map: Optional[Dict[str, str]] = None,
  gene_db_genes: Optional[Set[str]] = None,
) -> str:
  raw = str(gene).strip()
  if not raw:
    return raw

  cache_key = raw.lower()
  cached = _MAPPING_RESULT_CACHE.get(cache_key)
  if cached is not None:
    return cached

  if is_arabidopsis_tair_id(raw):
    result = _arabidopsis_tair_base(raw)
    _MAPPING_RESULT_CACHE[cache_key] = result
    return result

  sym_key = cache_key
  symbol_map = db8_map if db8_map is not None else _arabidopsis_db8_symbol_map(db8)
  hit = symbol_map.get(sym_key)
  if hit and is_arabidopsis_tair_id(hit):
    result = _arabidopsis_tair_base(hit)
    _MAPPING_RESULT_CACHE[cache_key] = result
    return result

  gene_db = gene_db_genes if gene_db_genes is not None else _arabidopsis_gene_db_genes()
  preferred = _pick_preferred_tair(_arabidopsis_tair_ids_from_ncbi(raw), gene_db)
  result = preferred if preferred else raw
  _MAPPING_RESULT_CACHE[cache_key] = result
  return result


def normalize_arabidopsis_marker_list(markers, *, db8=None) -> List[str]:
  if not markers:
    return []
  if all(is_arabidopsis_tair_id(marker) for marker in markers):
    return [_arabidopsis_tair_base(marker) for marker in markers]

  db8_map = _arabidopsis_db8_symbol_map(db8)
  needs_ncbi = any(
    not is_arabidopsis_tair_id(marker)
    and marker.strip().lower() not in db8_map
    and marker.strip().lower() not in _MAPPING_RESULT_CACHE
    for marker in markers
  )
  gene_db = _arabidopsis_gene_db_genes() if needs_ncbi else None
  return [
    map_arabidopsis_gene_to_tair(
      marker,
      db8_map=db8_map,
      gene_db_genes=gene_db,
    )
    for marker in markers
  ]


def normalize_arabidopsis_markers_str(markers: str, *, db8=None) -> str:
  if not markers or not str(markers).strip():
    return markers
  parts = [part.strip() for part in str(markers).split(",")]
  normalized = normalize_arabidopsis_marker_list(parts, db8=db8)
  return ", ".join(normalized)


def get_info_proteinAtlas(marker, db5, query):
  key = str(marker).lower().strip()
  mask = db5["Gene"].astype(str).str.lower().str.strip() == key
  record = db5[mask]
  
  doc = ""

  
  if query == "Gene description" and len(record["Gene description"]) != 0 and type(record["Gene description"].item()) == str:
    doc = f"{marker}: " + record["Gene description"].item()
  if query == "Protein class" and len(record["Protein class"]) != 0 and type(record["Protein class"].item()) == str:
    doc = f"{marker}'s protein class: " + record["Protein class"].item()
  if query == "Biological process" and len(record) > 0:
    val = record["Biological process"].iloc[0]
    if pd.notna(val):
      s = str(val).strip()
      if s:
        doc = f"{marker}'s biological process: " + s
  if query == "Molecular function" and len(record["Molecular function"]) != 0 and type(record["Molecular function"].item()) == str:
    doc = f"{marker}'s molecular function: " + record["Molecular function"].item()
  if query == "Tissue specificity" and len(record["RNA tissue specificity"]) != 0 and type(record["RNA tissue specificity"].item()) == str:
    doc += f"{marker}'s RNA tissue specificity: " + record["RNA tissue specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA tissue distribution"]) != 0 and type(record["RNA tissue distribution"].item()) == str:
    doc += f"{marker}'s RNA tissue distribution: " + record["RNA tissue distribution"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA single cell type specificity"]) != 0 and type(record["RNA single cell type specificity"].item()) == str:
    doc += f"{marker}'s RNA single cell type specificity: " + record["RNA single cell type specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA single cell type distribution"]) != 0 and type(record["RNA single cell type distribution"].item()) == str:
    doc += f"{marker}'s RNA single cell type distribution: " + record["RNA single cell type distribution"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA cancer specificity"]) != 0 and type(record["RNA cancer specificity"].item()) == str:
    doc += f"{marker}'s RNA cancer specificity: " + record["RNA cancer specificity"].item() + "\n"
  if query == "Tissue specificity" and len(record["RNA cancer distribution"]) != 0 and type(record["RNA cancer distribution"].item()) == str:
    doc += f"{marker}'s RNA cancer distribution: " + record["RNA cancer distribution"].item() + "\n"
  
  return doc


def get_uniprot_id(gene, species):
    url = "https://rest.uniprot.org/uniprotkb/search"
    
    if species == "mouse":
      species = "mus musculus"
    if species == "arabidopsis" or species == "Arabidopsis thaliana":
      species = "arabidopsis thaliana"
    
    tax_ids = {"human":9606, "mus musculus":10090, "arabidopsis thaliana":3702}
    
    tax_id = tax_ids[species.lower()]
    
    params = {
        "query": f"gene_exact:{gene} AND organism_id:{tax_id}",
        "fields": "accession,gene_primary,organism_name",
        "format": "json",
        "size": 1
    }
    
    try:
      r = requests.get(url, params=params)
      data = r.json()
    except Exception as e:
      return None

    if data.get("results"):
        result = data["results"][0]
        return result["primaryAccession"]
    else:
        return None



def scene_gene_desc_summary(state):
  markers = state["markers"].split(", ")
  gene_desc_summary_context = "[Gene's Description & Summary Part]\n"
  genetriever_flag = state["genetriever_flag"]
  species = state["species"]
  geneDB = get_gene_db_from_state(state)
  unknown_genes = state["unknown_genes"]

  if species.lower() == "arabidopsis":
    species = "Arabidopsis thaliana"
  db5 = state["db5"]
  genes_ids = {}
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      gene_desc_summary_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["Description"].item() + "\n"

      continue
    else:
      if genetriever_flag == False:
        unknown_genes.append(marker)
        continue
  
    gene_desc_summary_context += f"[Gene's Description & Summary]\n{marker}:\n"
    
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    
    params = {
       "db": "gene",
       "term": f"{marker}[gene] AND {species}[orgn]",
       "retmode": "json"
    }
    
    ncbi_desc = ""
    ncbi_summary = ""
    
    try:
      res = requests.get(esearch_url, params=params).json()
      gene_id = res["esearchresult"]["idlist"][0]
      genes_ids[marker] = gene_id
    
      efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    
      efetch = requests.get(efetch_url, params={"db":"gene", "id":gene_id, "retmode":"xml"})
    
      root = ElementTree.fromstring(efetch.content)
      
      if species == "Arabidopsis thaliana":
        desc = root.find(".//Prot-ref_desc")
        summary = root.find(".//Entrezgene_summary")
      else:
        desc = root.find(".//Gene-ref_desc")
        summary = root.find(".//Entrezgene_summary")

      try:
        ncbi_desc = desc.text
      except Exception as e:
        pass
    
      try:
        ncbi_summary = summary.text
      except Exception as e:
        pass
    except Exception as e:
      pass
    
    proteinAtlas_desc = ""
    proteinAtlas_class = ""
    if species == "human":  
      proteinAtlas_desc = get_info_proteinAtlas(marker, db5, "Gene description")
      proteinAtlas_class = get_info_proteinAtlas(marker, db5, "Protein class")
    
    try:
      uniprot_id = get_uniprot_id(marker, species)
      r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
      uniprot_docs = r.json()
    except Exception as e:
      uniprot_docs = None
      
    
    uniprot_desc = ""
    
    if uniprot_docs == None:
      pass
    else:
      try:
        protein_description = uniprot_docs['proteinDescription']
        if "recommendedName" in protein_description:
          uniprot_desc = f"{marker} encodes " + protein_description['recommendedName']['fullName']['value']

        else:
          uniprot_desc = f"{marker} encodes " + protein_description['submissionNames'][0]['fullName']['value']
    
      except Exception as e:
        pass
    
    gene_desc_summary_context += ncbi_desc + "\n" + ncbi_summary + "\n" + proteinAtlas_desc + "\n" + proteinAtlas_class + "\n" + uniprot_desc + "\n"

  pbar = state["pbar"]
  pbar.update(10)
  
  return {**state, "gene_desc_summary_context":gene_desc_summary_context, "genes_ids":genes_ids, "pbar":pbar, "unknown_genes":unknown_genes}
    
    
    
def scene_gene_function(state):
  markers = state["markers"].split(", ")
  gene_func_summary_context = "[Gene's function Part]\n"
  species = state["species"]
  genetriever_flag = state["genetriever_flag"]

  geneDB = get_gene_db_from_state(state)
  
  genes_ids = state["genes_ids"]
  db5 = state["db5"]
  llm = OpenAI()
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      gene_func_summary_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["Function"].item() + "\n"
      continue
    else:
      if genetriever_flag == False:
        continue

    gene_func_summary_context += f"[Gene's function]\n{marker}:\n"
    
    mygene_MF_text = ""
    
    try:
      gene_id = genes_ids[marker]
      fields = "go"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      go_text = mygene_response.text
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's biological process(BP), molecular function(MF), and cellular component(CC). You should extract MF contents only and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {go_text}

    """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-5.2",
        messages = messages,
        temperature = 0,
        )
    
      mygene_MF_text = completion.choices[0].message.content
    except Exception as e:
      pass
    
    proteinAtlas_mf = ""
    if species == "human":
      proteinAtlas_mf = get_info_proteinAtlas(marker, db5, "Molecular function")
    
    try:
      uniprot_id = get_uniprot_id(marker, species)
      r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
      uniprot_docs = r.json()
    except Exception as e:
      uniprot_id = None
      uniprot_docs = None
    
    uniprot_func = ""
    
    if uniprot_docs == None:
      pass
    else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "FUNCTION":
            for j in curPos['texts']:
              uniprot_func += f"{marker}'s function: " + j['value'] + "\n"
          elif curPos["commentType"] == "CATALYTIC ACTIVITY":
              uniprot_func += f"{marker}'s catalytic activity: " + curPos['reaction']['name'] + "\n"
              
      except Exception as e:
        pass
    
    gene_func_summary_context += mygene_MF_text + "\n" + proteinAtlas_mf + "\n" + uniprot_func + "\n"
  
  pbar = state["pbar"]
  pbar.update(10)
  return {**state, "gene_func_summary_context":gene_func_summary_context, "pbar":pbar}


def route_graph_construction(state):
  flag1 = state["graph_provided"]
  flag2 = state["genetriever_flag"]
  if flag1:
    return "Construct_Cell_Marker_Graph"
  # When genetriever is on, markers may be missing from geneDB but are NOT added to
  # unknown_genes (that list is only used when genetriever is off). scene_cell_relatedTo_gene
  # still needs marker_cell_graph for those markers — always build the graph first.
  if flag2:
    return "Construct_Cell_Marker_Graph"
  return "scene_cell_relatedTo_gene"



def scene_cell_relatedTo_gene(state):
  graph = state["marker_cell_graph"]
  markers = state["markers"].split(", ")
  gene_cell_summary_context = "[Cell types or Tissue related to Gene Part]\n"
  species = state["species"]
  genetriever_flag = state["genetriever_flag"]

  geneDB = get_gene_db_from_state(state)
  genes_ids = state["genes_ids"]
  db5 = state["db5"]
  
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      gene_cell_summary_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["Related cell type"].item() + "\n"
      continue
    else:
      if genetriever_flag == False:
        continue
      
    gene_cell_summary_context += f"[Cell types or Tissue related to Gene]\n{marker}:\n"
    
    cell_marker_text = f"{marker} was reported to be expressed in these kinds of cell types (There may be more cell types that {marker} can be expressed.): "

    # Isolated marker vertices are removed in construct_marker_graph; the name may be absent.
    if graph is None or marker not in graph.vs["name"]:
      neighbors_celltypes = []
    else:
      neighbors = graph.neighbors(marker)
      neighbors_celltypes = [graph.vs[idx]["name"] for idx in neighbors]

    cell_marker_text += " / ".join(neighbors_celltypes) if neighbors_celltypes else "(no cell types linked in marker graph after filtering)"
    
    proteinAtlas_tissue = ""
    
    if species == "human":
      proteinAtlas_tissue = get_info_proteinAtlas(marker, db5, "Tissue specificity")
    
    try:
      uniprot_id = get_uniprot_id(marker, species)
      r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
      uniprot_docs = r.json()
    except Exception as e:
      uniprot_id = None
      uniprot_docs = None
      
    uniprot_tissue = ""
    
    if uniprot_docs == None:
      pass
    else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "TISSUE SPECIFICITY":
            for j in curPos['texts']:
              uniprot_tissue += f"{marker}: " + j['value']
              
      except Exception as e:
        pass
    
    gene_cell_summary_context += cell_marker_text + "\n" + proteinAtlas_tissue + "\n" + uniprot_tissue + "\n"

  pbar = state["pbar"]
  pbar.update(10)
  return {**state, "gene_cell_summary_context":gene_cell_summary_context, "pbar":pbar}
  
  
def scene_pathway(state):
  gene_pathway_context = "[Pathway related to Gene Part]\n"
  llm = OpenAI()
  markers = state["markers"].split(", ")
  genes_ids = state["genes_ids"]
  db5 = state["db5"]
  species = state["species"]
  genetriever_flag = state["genetriever_flag"]

  geneDB = get_gene_db_from_state(state)
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      gene_pathway_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["Pathway"].item() + "\n"
      continue
    else:
      if genetriever_flag == False:
        continue
  
  
    if species == "arabidopsis":
      try:
        gene_pathway_context += f"[Pathway related to Gene]\n{marker}:\n"
        query=f"(gene:{marker}) AND (organism_id:3702)"
        url=f"https://rest.uniprot.org/uniprotkb/search?query={query}&fields=accession,gene_names,go_p"
        res = requests.get(url, headers={"Accept":"application/json"})
        data = res.json()
      
      except Exception as e:
        data = None
        
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene. You should extract biological process contents only and summarize it. Biological process is organized as the specific form (ex. value: P:cell division)
    
    [Gene]
    {marker}
    
    [Context]
    {data}"""
    
      messages = [{"role": "user", "content": prompt}]
      try:
        completion = llm.chat.completions.create(
          model = "gpt-5.2",
          messages = messages,
          temperature = 0,
        )
        pathway_text = completion.choices[0].message.content
        gene_pathway_context += pathway_text + "\n"
      except Exception as e:
        pathway_text = ""
      
      continue
      
      
    gene_pathway_context += f"[Pathway related to Gene]\n{marker}:\n"
    mygene_pathway_text = ""
    try:
      gene_id = genes_ids[marker]
    
      fields = "pathway"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      pathway_text = mygene_response.text
      
    
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's pathway. You should summarize it.

    [Gene]
    {marker}
    
    [Context]
    {pathway_text}
    
      """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-5.2",
        messages = messages,
        temperature = 0,
      )
    
      mygene_pathway_text = completion.choices[0].message.content
    
    except Exception as e:
      pass
    
    proteinAtlas_bp = ""
    if species == "human":
      proteinAtlas_bp = get_info_proteinAtlas(marker, db5, "Biological process")
    
    
    gene_pathway_context += mygene_pathway_text + "\n" + proteinAtlas_bp + "\n"

  pbar = state["pbar"]
  pbar.update(10)
  return {**state, "gene_pathway_context":gene_pathway_context, "pbar":pbar}

def scene_protein_location(state):
  gene_protein_loc_context = "[Gene's Protein Location Part]\n"
  llm = OpenAI()
  markers = state["markers"].split(", ")
  genes_ids = state["genes_ids"]
  species = state["species"]
  genetriever_flag = state["genetriever_flag"]

  geneDB = get_gene_db_from_state(state)
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      gene_protein_loc_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["Protein location"].item() + "\n"
      continue
    else:
      if genetriever_flag == False:
        continue
      
    gene_protein_loc_context += f"[Gene's Protein Location]\n{marker}:\n"
    mygene_cc_text = ""
    
    try:
      gene_id = genes_ids[marker]
    
      fields = "go"
      mygene_url = f"https://mygene.info/v3/gene/{gene_id}?fields={fields}"
      mygene_response = requests.get(mygene_url)
      go_text = mygene_response.text
    
      prompt = f"""You are an intelligent summarizer that summarizes the context of gene. You will be given the structured data about gene's biological process(BP), molecular function(MF), and cellular component(CC). You should extract CC contents only and summarize it.
    
    [Gene]
    {marker}
    
    [Context]
    {go_text}
    
      """
      messages = [{"role": "user", "content": prompt}]
    
      completion = llm.chat.completions.create(
        model = "gpt-5.2",
        messages = messages,
        temperature = 0,
      )
    
      mygene_cc_text = completion.choices[0].message.content
    
    except Exception as e:
      pass
    
    try:
      uniprot_id = get_uniprot_id(marker, species)
      r = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"Accept": "application/json"})
      uniprot_docs = r.json()
    except Exception as e:
      uniprot_id = None
      uniprot_docs = None
      
    uniprot_cc = ""
    
    if uniprot_docs == None:
      pass
    else:
      try:

        protein_comments = uniprot_docs['comments']
        
        for i in range(len(protein_comments)):
          curPos = protein_comments[i]
          
          if curPos["commentType"] == "SUBCELLULAR LOCATION":
            for j in curPos['subcellularLocations']:
              uniprot_cc += f"Proteins encoded by {marker} are located on " + j['location']['value']
              
      except Exception as e:
        pass
    
    
    
    gene_protein_loc_context += mygene_cc_text + "\n" + uniprot_cc + "\n"
  
  pbar = state["pbar"]
  pbar.update(10)
  
  
  return {**state, "gene_protein_loc_context":gene_protein_loc_context, "pbar":pbar}
  

def convert_id_to_name(ensemblID):
  mg = mygene.MyGeneInfo()
  buf_err = io.StringIO()
  with contextlib.redirect_stderr(buf_err):
    res = mg.query(ensemblID, scopes='ensembl.gene', fields='symbol')
  if len(res['hits']) != 0:
    try:
      return res['hits'][0]['symbol']
    except Exception as e:
      return ""
  else:
    return ""
    
    

def convert_name_to_id(name):
  mg = mygene.MyGeneInfo()
  buf_err = io.StringIO()
  with contextlib.redirect_stderr(buf_err):
    res = mg.query(name, scopes="symbol", fields="ensembl.gene", species="human")
  if len(res['hits']) != 0:
    try:
      return res['hits'][0]['ensembl']
    except Exception as e:
      return ""
  else:
    return ""
    
    

def p_p_interaction_info(state):
  ppi_context = "[Gene's Protein-Protein Interaction Part]\n"
  

  
  markers = state["markers"].split(", ")
  species = state["species"]
  db6 = state["db6"]
  genetriever_flag = state["genetriever_flag"]

  geneDB = get_gene_db_from_state(state)

  
  tax_ids = {"human":9606, "mouse":10090, "arabidopsis":3702}
  
  for marker in markers:
    if marker.lower() in geneDB["Gene"].str.lower().values:
      ppi_context += geneDB.loc[geneDB["Gene"].str.lower()==marker.lower()]["PP interaction"].item() + "\n"
      continue
    else:
      if genetriever_flag == False:
        continue
      
    ppi_context += f"[Gene's Protein-Protein Interaction]\n{marker}:\n"
    
    if species == "arabidopsis":
      pass
    else:
      marker_ids = convert_name_to_id(marker)
      if marker_ids == "":
        continue
      if type(marker_ids) == dict:
        marker_ids = [marker_ids]
      docs = []
      for marker_id in marker_ids:
        mask = (db6["ensembl_gene_id_1"] == marker_id['gene']) | (db6["ensembl_gene_id_2"] == marker_id['gene'])
        filter_db = db6[mask]
        partner_info = f"{marker} can interact with : "
  
        for record1, record2 in zip(filter_db['ensembl_gene_id_1'], filter_db['ensembl_gene_id_2']):
          if record1 == marker_id['gene']:
            partner_info += convert_id_to_name(record2)
            partner_info += " / "
          else:
            partner_info += convert_id_to_name(record1)
  
      ppi_context += partner_info + "\n"
    
    params = {
    "identifiers": marker,
    "species": tax_ids[species.lower()],
    "limit": 10,
    "required_score": 700,
    "format": "tsv"
    }
    
    try:
      res = requests.get("https://string-db.org/api/tsv/network", params=params)    
      llm = OpenAI()
      prompt = f"""You are an intelligent summarizer that summarizes the network of genes. You will be given table representing the interaction, so you should organize the network of them with natural sentences. You should only represent the gene as symbol. 
    
    [Interaction Tables]
    {res.text}
    
    """
    
      messages = [{"role": "user", "content": prompt}]
      completion = llm.chat.completions.create(
        model = "gpt-5.2",
         messages = messages,
         temperature = 0,
      )
      ppi_context += completion.choices[0].message.content + "\n"
    except Exception as e:
      completion = None
  
  pbar = state["pbar"]
  pbar.update(10)
  

  
  return {**state, "ppi_context":ppi_context, "pbar":pbar}

def route_unknown_genes(state):
  unknown_genes = state["unknown_genes"]
  if len(unknown_genes) == 0:
    return "select_candidates_from_llm"
  else:
    return "unknown_genes_info"
  
def unknown_genes_info(state):
  unknown_genes = state["unknown_genes"]
  unknown_genes_info = ""
  for gene in unknown_genes:
    unknown_genes_info += f"{gene} "
  
  print(unknown_genes_info + " are not found in the database. If you want to use these genes' information, please set genetriever_flag to True.")
  
  return {**state}
