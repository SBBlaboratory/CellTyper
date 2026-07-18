<h1 align="center"> CellTyper </h1> <br>
<p align="center">
  <a href="https://github.com/hyun-jin891/CellTyper">
    <img alt="CellTyper" title="CellTyper" src="https://i.postimg.cc/JnbjNcgV/celltyper-logo.jpg" width="450">
  </a>
</p>

<p align="center">
  LLM-based cell type annotation for single-cell marker genes<br>
  (v1: human · mouse · <i>Arabidopsis thaliana</i>)
</p>

<p align="center">
  <b>License:</b> Apache-2.0 (code) · see <a href="data_license.md">data_license.md</a> for third-party data
</p>

---

## Overview

**CellTyper** annotates scRNA-seq clusters from a marker gene table (Seurat `FindAllMarkers` or Scanpy `rank_genes_groups`).

For each cluster it:

1. Selects top marker genes  
2. Builds gene context (geneDB and optional live GeneTriever retrieval)  
3. Optionally builds a marker–cell graph  
4. Runs LLM self-consistency voting  
5. Writes an HTML report (prediction, contexts, consistency chart, graph)

---

## Requirements

- Python 3.10+
- An [OpenAI API key](https://platform.openai.com/api-keys)
- Reference data under `data/` (or auto-download via `CELLTYPER_DATA_BASE_URL` — see `celltyper_data.py`)
- Dependencies listed in `requirements.txt`

---

## Installation

```bash
git clone https://github.com/hyun-jin891/CellTyper.git
cd CellTyper

# recommended: create a dedicated environment
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

pip install -r requirements.txt
```

### `requirements.txt`

- Install **exactly** the pinned versions used for development/testing:
  ```bash
  pip install -r requirements.txt
  ```
- Prefer a **fresh virtual environment** (`venv` or `conda`) so versions do not clash with other projects.
- The file is a full dependency freeze (direct + transitive packages). Upgrading individual packages may break compatibility; if you change versions, re-test CellTyper end-to-end.
- After install, keep your OpenAI key available (see below) before calling `run_celltyper`.

---

## OpenAI API key

CellTyper needs an OpenAI key for annotation (and for GeneTriever LLM summaries when enabled).

**Option A — pass the key into `run_celltyper`:**

```python
from celltyper_public import run_celltyper

run_celltyper(
    markers="markers.csv",
    species="human",
    tissue="liver",
    condition="normal",
    openai_api_key="sk-...",
)
```

**Option B — set once per session:**

```python
from celltyper_public import set_openai_api_key, run_celltyper

set_openai_api_key("sk-...")
run_celltyper(
    markers="markers.csv",
    species="human",
    tissue="liver",
    condition="normal",
)
```

**Option C — environment variable:**

```bash
export OPENAI_API_KEY="sk-..."
```

```python
from celltyper_public import run_celltyper

run_celltyper(
    markers="markers.csv",
    species="human",
    tissue="liver",
    condition="normal",
)
```

If no key is available, `run_celltyper` raises a clear error.

---

## Quick start

```python
from celltyper_public import run_celltyper

result = run_celltyper(
    markers="path/to/markers.csv",
    species="human",          # human | mouse | arabidopsis
    tissue="liver",
    condition="normal",
    genetriever_flag=False,   # True: live retrieval for genes missing from geneDB
    graph_provided=True,      # include marker–cell graph in the workflow / report
    html_report_path="CellTyper_report.html",
    llm_model="gpt-5.2",
    max_workers=1,            # >1: annotate clusters in parallel
    openai_api_key="sk-...",  # or use set_openai_api_key / OPENAI_API_KEY
)

print(result["cluster_ids"])
print(result["celltypes"])
print(result["html_report_path"])
```

---

## Marker table formats

Gene symbols are preferred (Ensembl IDs like `ENSG…` / `ENSMUSG…` are dropped for human/mouse; *Arabidopsis* symbols are mapped to TAIR IDs automatically).

### Seurat / R (`FindAllMarkers`)

Required columns include:

| Column | Role |
|---|---|
| `gene` | gene symbol |
| `cluster` | cluster id |
| `p_val_adj` | adjusted p-value |
| `avg_log2FC` | log fold-change |
| `pct.1`, `pct.2` | detection fractions |

Filtering (summary): `p_val_adj < 0.05`, `pct.1 > 0.5`, rank by `avg_log2FC × (pct.1 − pct.2)`, then `avg_log2FC > 0.25`, top 10 per cluster.

### Scanpy (`rank_genes_groups`)

Export with `sc.get.rank_genes_groups_df` (use `pts=True` when possible). CellTyper accepts:

| Column | Notes |
|---|---|
| `names` or `gene` | gene id/symbol |
| `group` or `cluster` | cluster label |
| `pvals_adj` | adjusted p-value |
| `logfoldchanges` / `scores` | ranking |
| `pct_nz_group`, `pct_nz_reference` | optional; used like Seurat pct columns |

---

## Main parameters

| Parameter | Default | Description |
|---|---|---|
| `markers` | — | Path to marker CSV |
| `species` | — | `human`, `mouse`, or `arabidopsis` |
| `tissue` | — | Tissue name (must match supported list for that species) |
| `condition` | — | Free-text condition passed to the LLM |
| `genetriever_flag` | `False` | Live NCBI / MyGene / UniProt / STRING / … retrieval for geneDB misses |
| `graph_provided` | `False` | Build/use marker–cell graph (also controls report graph section) |
| `html_report_path` | `CellTyper_report.html` | Output HTML path (`None` to skip) |
| `llm_model` | `gpt-5.2` | OpenAI model for final annotation |
| `max_workers` | `1` | Max parallel clusters |
| `openai_api_key` | `None` | API key for this run (optional if env / `set_openai_api_key` already set) |

---

## Output

`run_celltyper` returns a dict:

```python
{
  "cluster_ids": [...],
  "celltypes": [...],           # lowercase keys used internally
  "html_report_path": "...",    # if HTML was written
}
```

The HTML report includes predicted labels, per-marker gene context, self-consistency (pie chart), LLM samples, and (when enabled) the marker–cell graph.

---

## Notes

- Final cell-type calls use **self-consistency** (multiple LLM samples + vote; ties are reported together).
- With `genetriever_flag=False`, context comes from bundled **geneDB** files when markers are present there.
- OpenAI prompt caching / Responses API paths are used for GPT-5-family models; other providers are not wired yet.

---

## License

**CellTyper source code** is released under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

Third-party **databases and reference data** have their own terms. See [`data_license.md`](data_license.md) for the list of resources (e.g. CellMarker, PanglaoDB, Human Protein Atlas, UniProt, STRING) and their licenses. Using or redistributing those data files remains subject to the original providers’ terms.

---

## Citation / related

Gene context retrieval builds on ideas from [GeneTriever](https://github.com/hyun-jin891/GeneTriever).
