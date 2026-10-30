#ST
"""
Agent tools: list_dir, read_file, inspect_h5ad, view_image, run_python, fetch_gene_map.
Each tool takes plain kwargs and returns a string (view_image returns image blocks);
errors are raised and turned into is_error tool results by the loop.
"""

# IMPORTS
import base64
import json
import re
import sys
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import anndata as ad

# src/01_preprocess starts with a digit, so it is not importable as a package
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "01_preprocess"))
from cs223_scRNAdataset_read_funs import PROJECT_ROOT, DATA_DIR
from gene_map import GENE_MAP_PATH, DOWNLOAD_NOTE, download_gene_map, load_gene_map

MAX_CHARS = 20000   # cap on any tool output so one call can't flood the context
MAX_LIST = 200
# Secrets and non-project dirs the agent must never see
BLOCKED = {".git", ".venv", ".idea", ".claude", "internal_docs", "config"}


def _resolve(path):
    """Project-relative path -> absolute Path, refusing anything outside the project or secret."""
    p = (PROJECT_ROOT / path).resolve()
    try:
        rel = p.relative_to(PROJECT_ROOT)
    except ValueError:
        raise ValueError(f"path is outside the project: {path}")
    if any(part in BLOCKED or part.startswith(".env") for part in rel.parts):
        raise ValueError(f"path is not accessible: {path}")
    return p


def _clip(text):
    if len(text) <= MAX_CHARS:
        return text
    return text[:MAX_CHARS] + f"\n... [truncated, {len(text) - MAX_CHARS} more chars]"


def list_dir(path="."):
    p = _resolve(path)
    if not p.is_dir():
        raise ValueError(f"not a directory: {path}")
    rows = []
    for child in sorted(p.iterdir()):
        if child.name in BLOCKED or child.name.startswith(".env"):
            continue
        rows.append(f"{child.name}/" if child.is_dir() else f"{child.name}  ({child.stat().st_size:,} bytes)")
    more = f"\n... [{len(rows) - MAX_LIST} more entries]" if len(rows) > MAX_LIST else ""
    return "\n".join(rows[:MAX_LIST]) + more if rows else "(empty)"


def read_file(path):
    p = _resolve(path)
    if not p.is_file():
        raise ValueError(f"not a file: {path}")
    try:
        return _clip(p.read_text())
    except UnicodeDecodeError:
        raise ValueError(f"{path} is binary; use inspect_h5ad for .h5ad files")


def _sample_metadata(p):
    """Row of data/metadata.json for Sample_<id>.h5ad, or None."""
    m = re.fullmatch(r"Sample_(\d+)", p.stem)
    meta_file = DATA_DIR / "metadata.json"
    if not m or not meta_file.exists():
        return None
    for row in json.loads(meta_file.read_text()):
        if str(row["sample_id"]) == m.group(1):
            return row
    return None


def inspect_h5ad(path, n_cells=2000):
    p = _resolve(path)
    if p.suffix != ".h5ad" or not p.is_file():
        raise ValueError(f"not an .h5ad file: {path}")

    # backed mode: read structure without loading the whole matrix into memory
    adata = ad.read_h5ad(p, backed="r")
    try:
        out = [f"file: {path}", f"shape: {adata.n_obs} cells x {adata.n_vars} genes"]
        meta = _sample_metadata(p)
        if meta:
            out.append(f"metadata.json: age={meta['age']}, tissue={meta['tissue']}, sex={meta['sex']}, "
                       f"disease={meta['disease_status']}, format={meta['format']}")

        out.append(f"obs columns: {list(adata.obs.columns)}")
        out.append(f"var columns: {list(adata.var.columns)}")
        out.append(f"obs_names: {list(adata.obs_names[:3])} ...")
        out.append(f"var_names: {list(adata.var_names[:5])} ...")
        # layers.keys() is unreliable in backed mode, so read the names from the file
        layer_names = list(adata.file._file["layers"].keys()) if "layers" in adata.file._file else []
        out.append(f"layers: {layer_names}; obsm: {list(adata.obsm.keys())}; "
                   f"uns: {list(adata.uns.keys())}")

        # Look at the first n_cells only; enough to tell raw counts from normalized data
        k = min(n_cells, adata.n_obs)
        X = adata.X[:k]
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(np.asarray(X))
        vals = X.data
        is_int = bool(np.all(vals == np.round(vals))) if vals.size else True
        counts = np.asarray(X.sum(axis=1)).ravel()
        genes = np.diff(X.indptr)
        # var_names can be Ensembl IDs; MT- genes are only visible in the symbol column
        symbols = adata.var["gene_symbols"] if "gene_symbols" in adata.var else adata.var_names
        mt = np.array([str(n).upper().startswith("MT-") for n in symbols])
        if GENE_MAP_PATH.exists():   # preferred: chromosome from the gene map
            table = load_gene_map()
            mt = np.asarray(adata.var_names.isin(table.index[table["chrom"] == "MT"]))
        mt_frac = np.asarray(X[:, mt].sum(axis=1)).ravel() / np.maximum(counts, 1)
        q = lambda a: " / ".join(f"{v:.0f}" for v in np.percentile(a, [5, 50, 95]))

        out.append(f"--- X stats over first {k} cells ---")
        # backed sparse X is not a scipy matrix, but it exposes .format ('csr'/'csc')
        out.append(f"X: {getattr(adata.X, 'format', 'dense')} {adata.X.dtype}; "
                   f"min={vals.min() if vals.size else 0}, max={vals.max() if vals.size else 0}; "
                   f"integer-valued={is_int} (True suggests raw counts)")
        out.append(f"total counts per cell (p5/p50/p95): {q(counts)}")
        out.append(f"genes detected per cell (p5/p50/p95): {q(genes)}")
        if mt.sum() == 0:
            out.append("mito fraction: n/a (genes are Ensembl IDs without symbols and no gene map found; "
                       "call fetch_gene_map, then add_gene_info in run_python)")
        else:
            out.append("mito fraction (p5/p50/p95): "
                       + " / ".join(f"{v:.3f}" for v in np.percentile(mt_frac, [5, 50, 95]))
                       + f"  ({mt.sum()} MT- genes)")
        return "\n".join(out)
    finally:
        adata.file.close()


# Tool definitions sent to the API. eager_input_streaming: stream inputs as they are generated.
TOOLS = [
    {
        "name": "list_dir",
        "description": "List files and folders in a project directory (path relative to the project root).",
        "input_schema": {
            "type": "object",
            "properties": {"path": {"type": "string", "description": "e.g. 'data'. Default '.'"}},
        },
        "eager_input_streaming": True,
    },
    {
        "name": "read_file",
        "description": "Read a UTF-8 text file in the project (docs, metadata.json, scripts). Output is truncated at ~20k chars.",
        "input_schema": {
            "type": "object",
            "properties": {"path": {"type": "string", "description": "e.g. 'data/metadata.json'"}},
            "required": ["path"],
        },
        "eager_input_streaming": True,
    },
    {
        "name": "inspect_h5ad",
        "description": (
            "Summarize an .h5ad file without loading it all: shape, obs/var columns, layers, whether X looks "
            "like raw counts, QC-style stats (counts, genes, mito fraction) and the sample's metadata.json row."
        ),
        "input_schema": {
            "type": "object",
            "properties": {"path": {"type": "string", "description": "e.g. 'data/Sample_3444.h5ad'"}},
            "required": ["path"],
        },
        "eager_input_streaming": True,
    },
]

MEDIA_TYPES = {".png": "image/png", ".jpg": "image/jpeg", ".jpeg": "image/jpeg"}
MAX_IMAGE_BYTES = 5_000_000   # API image size limit


def view_image(path):
    """Returns the image as a content block list so Claude actually sees it."""
    p = _resolve(path)
    if p.suffix.lower() not in MEDIA_TYPES or not p.is_file():
        raise ValueError(f"not a .png/.jpg file: {path}")
    if p.stat().st_size > MAX_IMAGE_BYTES:
        raise ValueError(f"image is {p.stat().st_size:,} bytes, over the 5 MB limit; save it smaller (lower dpi)")
    data = base64.standard_b64encode(p.read_bytes()).decode()
    return [{"type": "image", "source": {"type": "base64", "media_type": MEDIA_TYPES[p.suffix.lower()], "data": data}}]


def fetch_gene_map(ask):
    """Builds data/gene_map.tsv if missing; downloads only after the user says yes."""
    if GENE_MAP_PATH.exists():
        return f"{GENE_MAP_PATH.name} already exists. Use add_gene_info(adata) in run_python."
    if not ask(f"fetch_gene_map wants to download: {DOWNLOAD_NOTE}"):
        raise RuntimeError("the user declined the download")
    table = download_gene_map()
    return f"Wrote {GENE_MAP_PATH.name} ({len(table)} genes). Use add_gene_info(adata) in run_python."


TOOLS += [
    {
        "name": "fetch_gene_map",
        "description": (
            "Create the Ensembl ID -> gene symbol/chromosome table (data/gene_map.tsv) if it is missing. "
            "The data only has ENSG IDs, so call this before you need symbols or mitochondrial genes. "
            "The user is asked to approve the download. If the table already exists it does nothing."
        ),
        "input_schema": {"type": "object", "properties": {}},
    },
    {
        "name": "run_python",
        "description": (
            "Run Python code in a persistent kernel (working dir = project root, matplotlib backend = Agg). "
            "Variables and imports survive between calls, so load data once. Only stdout, the last expression "
            "and errors come back, so print what you need to see. To look at a plot, save it with plt.savefig "
            "(a .png under output/ or analysis_results/) and then call view_image."
        ),
        "input_schema": {
            "type": "object",
            "properties": {
                "code": {"type": "string"},
                "timeout": {"type": "integer", "description": "seconds, default 120, max 600"},
            },
            "required": ["code"],
        },
        "eager_input_streaming": True,
    },
    {
        "name": "view_image",
        "description": "Look at a .png or .jpg file in the project (e.g. a plot you saved). Max 5 MB.",
        "input_schema": {
            "type": "object",
            "properties": {"path": {"type": "string", "description": "e.g. 'output/umap.png'"}},
            "required": ["path"],
        },
        "eager_input_streaming": True,
    },
]

_FUNCS = {"list_dir": list_dir, "read_file": read_file, "inspect_h5ad": inspect_h5ad, "view_image": view_image}


def run_tool(name, tool_input, kernel=None, ask=None):
    """Dispatch one tool call; raises on bad input so the loop can report it as is_error.
    `ask(question) -> bool` is how tools get the user's y/N."""
    if name == "fetch_gene_map":
        return fetch_gene_map(ask)
    if name == "run_python":
        timeout = min(int(tool_input.get("timeout", 120)), 600)
        return kernel.execute(tool_input["code"], timeout=timeout)
    if name not in _FUNCS:
        raise ValueError(f"unknown tool: {name}")
    result = _FUNCS[name](**tool_input)
    return _clip(result) if isinstance(result, str) else result
