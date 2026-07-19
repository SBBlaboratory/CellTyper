<h1 align="center"> CellTyperv1 </h1> <br>
<p align="center">
  <a href="https://github.com/hyun-jin891/CellTyper">
    <img alt="CellTyper" title="CellTyper" src="https://i.postimg.cc/JnbjNcgV/celltyper-logo.jpg" width="450">
  </a>
</p>

<p align="center">
  Cell-Type Annotation LLM Agent for single-cell marker genes<br>
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
2. Builds gene context (GeneDB and optional live GeneTriever retrieval)  
3. Optionally builds a marker–cell graph  
4. Runs LLM self-consistency voting  
5. Writes an HTML report (prediction, contexts, consistency chart, graph)

> **We strongly recommend reviewing the HTML report** (`CellTyper_report.html` by default) as the primary way to interpret CellTyper results.  
> The returned `celltypes` list alone omits gene context, self-consistency votes, LLM samples, and the marker–cell graph that explain each prediction.

---

## Requirements

- Python 3.10+
- An [OpenAI API key](https://platform.openai.com/api-keys)
- Reference data under `data/` (or auto-download via `CELLTYPER_DATA_BASE_URL` — see `celltyper_data.py`)
- Dependencies listed in `requirements.txt`

---

## Installation

```bash
git clone https://github.com/SBBlaboratory/CellTyper.git
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
    species="human",          # human | mouse | Arabidopsis
    tissue="liver",
    condition="normal",
    genetriever_flag=False,   # True: live retrieval for genes missing from GeneDB
    graph_provided=True,      # include marker–cell graph in the workflow / report
    html_report_path="CellTyper_report.html",
    llm_model="gpt-5.2",
    max_workers=3,            # recommended >1 to cut wall-clock time (parallel clusters)
    openai_api_key="sk-...",  # or use set_openai_api_key / OPENAI_API_KEY
)

print(result["cluster_ids"])
print(result["celltypes"])
print(result["html_report_path"])  # open this HTML report to review full results
```

**Important:** Treat the HTML report as the main result. Console / dict outputs are a short summary only.

**Performance tip:** Set `max_workers` **greater than 1** (e.g. `2`–`4`) so multiple clusters annotate in parallel and **reduce total computation / wall-clock time**. The default remains `1` for safer sequential runs; raise it when you have several clusters and acceptable OpenAI rate limits.

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

Filtering (summary):

- `pvals_adj < 0.05` (when the column is present)
- if `pct_nz_group` is present: `pct_nz_group > 0.5`
- ranking:
  - with `pct_nz_group` + `pct_nz_reference` + `logfoldchanges`:  
    `logfoldchanges × (pct_nz_group − pct_nz_reference)` (descending)
  - else prefer `scores`, then `logfoldchanges`
- per cluster: `logfoldchanges > 0.25` when available (otherwise score-ranked rows), then **top 10**

---

## Supported species and tissues

`species` and `tissue` must match these names **exactly** (shown in lowercase; `run_celltyper` lowercases inputs).

### `human`

`adipose tissue`, `adrenal gland`, `axilla`, `blood`, `bladder organ`, `bone marrow`, `brain`, `breast`, `colon`, `embryo`, `endocrine gland`, `esophagus`, `exocrine gland`, `eye`, `fallopian tube`, `heart`, `intestine`, `kidney`, `lamina propria`, `large intestine`, `liver`, `lung`, `lymph node`, `mucosa`, `musculature`, `nose`, `omentum`, `ovary`, `pancreas`, `pleural fluid`, `placenta`, `prostate gland`, `respiratory system`, `saliva`, `skeletal system`, `skin of body`, `small intestine`, `spleen`, `stomach`, `tongue`, `urinary bladder`, `uterus`, `vasculature`

### `mouse`

`kidney`, `adipose tissue`, `blood`, `bone marrow`, `brain`, `colon`, `embryo`, `endocrine gland`, `exocrine gland`, `eye`, `heart`, `large intestine`, `liver`, `lung`, `lymph node`, `mucosa`, `musculature`, `ovary`, `pancreas`, `prostate gland`, `respiratory system`, `skeletal system`, `skin of body`, `small intestine`, `spleen`, `tongue`, `urethra`, `urinary bladder`, `vasculature`

### `arabidopsis`

`shoot`, `seedling`, `root`, `leaf`, `fruit`, `flower`, `shoot apex`, `stem`

---

## Main parameters

| Parameter | Default | Description |
|---|---|---|
| `markers` | — | Path to marker CSV |
| `species` | — | `human`, `mouse`, or `arabidopsis` |
| `tissue` | — | Tissue name (must match the list above for that species) |
| `condition` | — | Free-text condition passed to the LLM |
| `genetriever_flag` | `False` | Live NCBI / MyGene / UniProt / STRING / … retrieval for GeneDB misses |
| `graph_provided` | `True` | Build/use marker–cell graph (also controls report graph section) |
| `html_report_path` | `CellTyper_report.html` | Output HTML path (`None` to skip) |
| `llm_model` | `gpt-5.2` | OpenAI model for final annotation |
| `max_workers` | `1` | Max parallel clusters. **Recommended `>1`** (e.g. 2–4) to save wall-clock time when annotating many clusters; watch API rate limits |
| `openai_api_key` | `None` | API key for this run (optional if env / `set_openai_api_key` already set) |

---

## Output

> **Recommended:** Always open and inspect the HTML report after a run.  
> It is the fullest view of CellTyper’s decisions and supporting evidence.

`run_celltyper` returns a dict:

```python
{
  "cluster_ids": [...],
  "celltypes": [...],           # lowercase keys used internally
  "html_report_path": "...",    # if HTML was written
}
```

The HTML report includes predicted labels, per-marker gene context, self-consistency (pie chart), LLM samples, and (when enabled) the marker–cell graph. Prefer this report over relying only on `celltypes` in the return value.

### Report preview

Example HTML report layout (scroll through each section in the generated file):

<p align="center">
  <img src="https://i.postimg.cc/63mBrXt1/seukeulinsyas-2026-07-19-112629.png" alt="CellTyper report — cluster header, markers, and gene context" width="720">
</p>
<p align="center"><em>Cluster summary, input markers, and expandable gene context</em></p>

<p align="center">
  <img src="https://i.postimg.cc/4dSX6gsM/seukeulinsyas-2026-07-19-112646.png" alt="CellTyper report — self-consistency and LLM responses" width="720">
</p>
<p align="center"><em>Self-consistency pie chart and LLM sample responses</em></p>

<p align="center">
  <img src="https://i.postimg.cc/d1xJ2Ywp/seukeulinsyas-2026-07-19-112700.png" alt="CellTyper report — marker–cell graph" width="720">
</p>
<p align="center"><em>Marker–cell graph (when <code>graph_provided=True</code>)</em></p>

---

## Notes

- Final cell-type calls use **self-consistency** (multiple LLM samples + vote; ties are reported together).
- With `genetriever_flag=False`, context comes from bundled **GeneDB** files when markers are present there.
- **Prefer `max_workers > 1`** for multi-cluster runs to shorten overall runtime (clusters processed in parallel). Start modestly (e.g. 2–3) if you hit rate limits.
- **OpenAI prompt caching / Responses API paths are used for GPT family models; other providers are not wired yet.**

---

## License

**CellTyper source code** is released under the [Apache License 2.0](https://www.apache.org/licenses/LICENSE-2.0).

Third-party **databases and reference data** have their own terms. See [`data_license.md`](data_license.md) for the list of resources (e.g. CellMarker, PanglaoDB, Human Protein Atlas, UniProt, STRING) and their licenses. Using or redistributing those data files remains subject to the original providers’ terms.

---

## Citation / related

Gene context retrieval builds on ideas from [GeneTriever](https://github.com/hyun-jin891/GeneTriever).


## Authors

- **Hyunjin Cho** — Lead developer, methodology, software implementation  
  [hyun-jin891](https://github.com/hyun-jin891)

- **Hyobin Jeong** — Supervision, project administration  
  [jeongdo801](https://github.com/jeongdo801)
