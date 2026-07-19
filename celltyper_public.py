"""Public API for programmatic CellTyper runs and HTML reports."""
from celltyper_data import data_file, ensure_celltyper_data

ensure_celltyper_data()

import base64
import os
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import pandas as pd
from tqdm import tqdm

import CellTyper
import pronto

from cell_marker_graph import construct_marker_graph
from celltyper_report import (
    build_html_report,
    build_per_marker_gene_contexts,
    marker_cell_graph_image_payload,
    write_html_report,
)
from geneRetriever import normalize_arabidopsis_marker_list

## Your Dataset ##
db1 = pd.read_excel(data_file("Cell_marker_Seq.xlsx"))
db2 = pd.read_csv(data_file("PanglaoDB_markers_27_Mar_2020.tsv"), sep="\t")
db3 = pd.read_csv(data_file("singleCellBase_dataset.txt"), sep="\t")
db5 = pd.read_csv(data_file("proteinatlas.tsv"), sep="\t")
db6 = pd.read_csv(data_file("interaction_consensus.tsv"), sep="\t")
ont = pronto.Ontology(data_file("cl-basic.obo"))
db8 = pd.read_excel(data_file("arabidopsis_thaliana.marker_fd.xlsx"))
## Your Dataset ##

_print_lock = threading.Lock()


def set_openai_api_key(api_key: str) -> None:
    """Set the OpenAI API key used by CellTyper / GeneTriever LLM calls.

    Equivalent to exporting ``OPENAI_API_KEY`` in the environment. Prefer
    passing ``openai_api_key=...`` to :func:`run_celltyper`, or calling this
    once at the start of a session.
    """
    key = "" if api_key is None else str(api_key).strip()
    if not key:
        raise ValueError("openai_api_key must be a non-empty string")
    os.environ["OPENAI_API_KEY"] = key


def _ensure_openai_api_key(openai_api_key=None) -> None:
    if openai_api_key is not None:
        set_openai_api_key(openai_api_key)
        return
    if not (os.environ.get("OPENAI_API_KEY") or "").strip():
        raise ValueError(
            "OpenAI API key is not set. Pass openai_api_key=... to run_celltyper(), "
            "call set_openai_api_key(...), or set the OPENAI_API_KEY environment variable."
        )


def _normalize_scanpy_marker_columns(markerTable: pd.DataFrame) -> pd.DataFrame:
    """Map Scanpy rank_genes_groups_df column names to CellTyper's gene/cluster schema."""
    out = markerTable.copy()
    if "gene" not in out.columns and "names" in out.columns:
        out["gene"] = out["names"]
    if "cluster" not in out.columns and "group" in out.columns:
        out["cluster"] = out["group"]
    return out


def _prepare_marker_table(markerTable: pd.DataFrame):
    """Filter/rank markers for Seurat (R) or Scanpy tables. Returns (table, fc_col, threshold)."""
    # Seurat / FindAllMarkers
    if "p_val_adj" in markerTable.columns:
        table = markerTable.copy()
        table = table[table["p_val_adj"] < 0.05]
        table = table[table["pct.1"] > 0.5]
        table["pct.1-pct.2"] = table["pct.1"] - table["pct.2"]
        table["log2FC * pct_diff"] = table["avg_log2FC"] * table["pct.1-pct.2"]
        table = table.sort_values("log2FC * pct_diff", ascending=False)
        return table, "avg_log2FC", 0.25

    # Scanpy / rank_genes_groups (+ optional pts=True proportions)
    scanpy_markers = (
        "pvals_adj" in markerTable.columns
        or "logfoldchanges" in markerTable.columns
        or "scores" in markerTable.columns
    )
    if scanpy_markers:
        table = _normalize_scanpy_marker_columns(markerTable)
        if "gene" not in table.columns:
            raise ValueError(
                "Scanpy marker table needs a gene column "
                "(use 'names' from rank_genes_groups_df or rename to 'gene')"
            )
        if "cluster" not in table.columns:
            raise ValueError(
                "Scanpy marker table needs a cluster column "
                "(use 'group' from rank_genes_groups_df or rename to 'cluster')"
            )
        if "pvals_adj" in table.columns:
            table = table[table["pvals_adj"] < 0.05]
        if "pct_nz_group" in table.columns:
            table = table[table["pct_nz_group"] > 0.5]
            if "logfoldchanges" in table.columns and "pct_nz_reference" in table.columns:
                table["pct.1-pct.2"] = table["pct_nz_group"] - table["pct_nz_reference"]
                table["log2FC * pct_diff"] = (
                    table["logfoldchanges"] * table["pct.1-pct.2"]
                )
                table = table.sort_values("log2FC * pct_diff", ascending=False)
            elif "logfoldchanges" in table.columns:
                table = table.sort_values("logfoldchanges", ascending=False)
            elif "scores" in table.columns:
                table = table.sort_values("scores", ascending=False)
        elif "scores" in table.columns:
            table = table.sort_values("scores", ascending=False)
        elif "logfoldchanges" in table.columns:
            table = table.sort_values("logfoldchanges", ascending=False)
        else:
            raise ValueError(
                "Scanpy marker table needs 'scores' and/or 'logfoldchanges' for ranking"
            )

        if "logfoldchanges" in table.columns:
            return table, "logfoldchanges", 0.25
        # logreg-style export has scores only; keep top-ranked rows after sort
        return table, "scores", float("-inf")

    raise ValueError(
        "Unrecognized marker table. Expected Seurat columns "
        "(p_val_adj, avg_log2FC, pct.1, pct.2) or Scanpy columns "
        "(pvals_adj, logfoldchanges/scores; optional pct_nz_group, pct_nz_reference)"
    )


def _build_cluster_jobs(markerTable, fc_col, threshold, species):
    jobs = []
    for index, val in enumerate(markerTable["cluster"].unique()):
        if str(val).isdigit():
            df = markerTable[markerTable["cluster"] == int(val)]
            cluster_label = str(val)
            is_named_cluster = False
        else:
            df = markerTable[markerTable["cluster"] == val]
            cluster_label = str(val)
            is_named_cluster = True

        df = df[df[fc_col] > threshold]
        if species == "mouse":
            df = df[~df["gene"].astype(str).str.startswith("ENSMUSG")]
        df = df.head(10)
        marker_gene_list = [str(x) for x in df["gene"].tolist()]
        if species == "arabidopsis":
            marker_gene_list = normalize_arabidopsis_marker_list(
                marker_gene_list, db8=db8
            )
        markers_str = ", ".join(marker_gene_list)

        jobs.append(
            {
                "index": index,
                "val": val,
                "cluster_label": cluster_label,
                "is_named_cluster": is_named_cluster,
                "marker_gene_list": marker_gene_list,
                "markers_str": markers_str,
            }
        )
    return jobs


def _process_cluster(
    job,
    *,
    species,
    tissue,
    condition,
    genetriever_flag,
    graph_provided,
    llm_model,
    progress_disable,
    announce_start,
):
    cluster_label = job["cluster_label"]
    marker_gene_list = job["marker_gene_list"]
    markers_str = job["markers_str"]

    if announce_start:
        with _print_lock:
            print(f"Processing cluster {cluster_label}")
            print()

    fin_states = CellTyper.app.invoke(
        {
            "db1": db1,
            "db2": db2,
            "db3": db3.copy(),
            "db5": db5,
            "db6": db6,
            "db8": db8,
            "ont": ont,
            "marker_cell_graph": None,
            "celltype_candidates": [],
            "fin_celltype": None,
            "pbar": None,
            "cluster_id": cluster_label,
            "progress_disable": progress_disable,
            "markers": markers_str,
            "species": species,
            "tissue": tissue,
            "condition": condition,
            "gene_desc_summary_context": "",
            "genes_ids": {},
            "gene_func_summary_context": "",
            "gene_cell_summary_context": "",
            "gene_pathway_context": "",
            "gene_protein_loc_context": "",
            "ppi_context": "",
            "genetriever_flag": genetriever_flag,
            "unknown_genes": [],
            "llm_model": "" if llm_model is None else str(llm_model).strip(),
            "graph_provided": graph_provided,
        },
        config={"recursion_limit": 10000},
    )

    res_celltype = fin_states["fin_celltype"]
    res_display = fin_states.get("fin_celltype_display") or res_celltype

    graph_b64 = None
    graph_mime = None
    graph_backend = None
    if graph_provided:
        g = fin_states.get("marker_cell_graph")
        if g is None:
            try:
                g = construct_marker_graph(fin_states)["marker_cell_graph"]
            except Exception:
                try:
                    st = {
                        **fin_states,
                        "pbar": tqdm(total=100, disable=True, leave=False),
                    }
                    g = construct_marker_graph(st)["marker_cell_graph"]
                except Exception:
                    g = None
        if g is not None:
            payload = marker_cell_graph_image_payload(g)
            if payload:
                raw, graph_mime, graph_backend = payload
                graph_b64 = base64.b64encode(raw).decode("ascii")

    cluster_row = {
        "cluster_id": cluster_label,
        "marker_genes": marker_gene_list,
        "predicted": res_celltype,
        "predicted_display": res_display,
        "llm_sample_texts": fin_states.get("llm_sample_texts") or [],
        "llm_parsed_celltypes": fin_states.get("llm_parsed_celltypes") or [],
        "llm_consistency_counts": fin_states.get("llm_consistency_counts") or {},
        "llm_consistency_scores": fin_states.get("llm_consistency_scores") or {},
        "gene_desc_summary_context": fin_states.get("gene_desc_summary_context") or "",
        "gene_func_summary_context": fin_states.get("gene_func_summary_context") or "",
        "gene_cell_summary_context": fin_states.get("gene_cell_summary_context") or "",
        "gene_pathway_context": fin_states.get("gene_pathway_context") or "",
        "gene_protein_loc_context": fin_states.get("gene_protein_loc_context") or "",
        "ppi_context": fin_states.get("ppi_context") or "",
        "marker_gene_contexts": build_per_marker_gene_contexts(
            marker_gene_list,
            species=species,
            gene_desc_summary_context=fin_states.get("gene_desc_summary_context") or "",
            gene_func_summary_context=fin_states.get("gene_func_summary_context") or "",
            gene_cell_summary_context=fin_states.get("gene_cell_summary_context") or "",
            gene_pathway_context=fin_states.get("gene_pathway_context") or "",
            gene_protein_loc_context=fin_states.get("gene_protein_loc_context") or "",
            ppi_context=fin_states.get("ppi_context") or "",
        ),
        "graph_png_b64": graph_b64,
        "graph_image_mime": graph_mime,
        "graph_export_backend": graph_backend,
    }

    return {
        "index": job["index"],
        "val": job["val"],
        "cluster_label": cluster_label,
        "is_named_cluster": job["is_named_cluster"],
        "res_celltype": res_celltype,
        "res_display": res_display,
        "cluster_row": cluster_row,
    }


def run_celltyper(
    markers,
    species,
    tissue,
    condition,
    genetriever_flag=False,
    graph_provided=True,
    html_report_path="CellTyper_report.html",
    llm_model="gpt-5.2",
    max_workers=1,
    openai_api_key=None,
):
    _ensure_openai_api_key(openai_api_key)

    markerTable = pd.read_csv(markers)

    species = species.lower()
    tissue = tissue.lower()
    condition = condition.lower()

    if max_workers < 1:
        raise ValueError("max_workers must be >= 1")

    if species not in ["human", "mouse", "arabidopsis"]:
        raise ValueError("Invalid species. Please choose from: human, mouse, arabidopsis")

    if species == "human":
        if tissue not in ["adipose tissue", "adrenal gland", "axilla", "blood", "bladder organ", "bone marrow", "brain", "breast", "colon", "embryo", "endocrine gland", "esophagus", "exocrine gland", "eye", "fallopian tube", "heart", "intestine", "kidney", "lamina propria", "large intestine", "liver", "lung", "lymph node", "mucosa", "musculature", "nose", "omentum", "ovary", "pancreas", "pleural fluid", "placenta", "prostate gland", "respiratory system", "saliva", "skeletal system", "skin of body", "small intestine", "spleen", "stomach", "tongue", "urinary bladder", "uterus", "vasculature"]:
            raise ValueError("Invalid tissue for human. Please choose from: adipose tissue, adrenal gland, axilla, blood, bladder organ, bone marrow, brain, breast, colon, embryo, endocrine gland, esophagus, exocrine gland, eye, fallopian tube, heart, intestine, kidney, lamina propria, large intestine, liver, lung, lymph node, mucosa, musculature, nose, omentum, ovary, pancreas, pleural fluid, placenta, prostate gland, respiratory system, saliva, skeletal system, skin of body, small intestine, spleen, stomach, tongue, urinary bladder, uterus, vasculature")
    elif species == "mouse":
        if tissue not in ["kidney", "adipose tissue", "blood", "bone marrow", "brain", "colon", "embryo", "endocrine gland", "exocrine gland", "eye", "heart", "large intestine", "liver", "lung", "lymph node", "mucosa", "musculature", "ovary", "pancreas", "prostate gland", "respiratory system", "skeletal system", "skin of body", "small intestine", "spleen", "tongue", "urethra", "urinary bladder", "vasculature"]:
            raise ValueError("Invalid tissue for mouse. Please choose from: kidney, adipose tissue, blood, bone marrow, brain, colon, embryo, endocrine gland, exocrine gland, eye, heart, large intestine, liver, lung, lymph node, mucosa, musculature, ovary, pancreas, prostate gland, respiratory system, skeletal system, skin of body, small intestine, spleen, tongue, urethra, urinary bladder, vasculature")
    elif species == "arabidopsis":
        if tissue not in ["shoot", "seedling", "root", "leaf", "fruit", "flower", "shoot apex", "stem"]:
            raise ValueError("Invalid tissue for arabidopsis. Please choose from: shoot, seedling, root, leaf, fruit, flower, shoot apex, stem")

    markerTable, fc_col, threshold = _prepare_marker_table(markerTable)

    jobs = _build_cluster_jobs(markerTable, fc_col, threshold, species)
    if not jobs:
        raise ValueError("No clusters found in markerTable after filtering")

    n_clusters = len(jobs)
    workers = min(max_workers, n_clusters)
    parallel = workers > 1
    progress_disable = parallel

    print(
        f"CellTyper: {n_clusters} cluster(s), "
        f"max_workers={max_workers} (using {workers} parallel slot(s))."
    )
    if parallel:
        print("Per-cluster progress bars are disabled during parallel runs.")
    print()

    shared_kw = dict(
        species=species,
        tissue=tissue,
        condition=condition,
        genetriever_flag=genetriever_flag,
        graph_provided=graph_provided,
        llm_model=llm_model,
        progress_disable=progress_disable,
        announce_start=not parallel,
    )

    results = []
    if parallel:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            future_to_job = {
                executor.submit(_process_cluster, job, **shared_kw): job for job in jobs
            }
            with tqdm(
                total=n_clusters,
                desc="Clusters completed",
                unit="cluster",
                bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt}",
            ) as overall:
                for future in as_completed(future_to_job):
                    result = future.result()
                    results.append(result)
                    overall.update(1)
                    tqdm.write(
                        f"Cluster {result['cluster_label']}: {result['res_display']}"
                    )
    else:
        for job in jobs:
            results.append(_process_cluster(job, **shared_kw))

    results.sort(key=lambda item: item["index"])

    result_celltype = [item["res_celltype"] for item in results]
    cluster_rows = [item["cluster_row"] for item in results]

    flag = any(item["is_named_cluster"] for item in results)
    cids = [item["cluster_label"] for item in results if item["is_named_cluster"]]
    vals = [item["val"] for item in results if not item["is_named_cluster"]]

    print()
    print("Summary")
    print("-------")
    for item in results:
        print(f"Cluster {item['cluster_label']}: {item['res_display']}")

    cluster_ids = cids if flag else vals

    out = {
        "cluster_ids": cluster_ids,
        "celltypes": result_celltype,
    }
    if html_report_path:
        html_doc = build_html_report(
            markers_path=str(markers),
            species=species,
            tissue=tissue,
            condition=condition,
            clusters=cluster_rows,
        )
        write_html_report(html_report_path, html_doc)
        out["html_report_path"] = html_report_path

    return out
