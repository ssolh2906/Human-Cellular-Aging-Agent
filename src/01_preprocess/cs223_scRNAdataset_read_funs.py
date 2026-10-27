#NAME: cs223_scRNAdataset_read_funs.py

# IMPORTS
import os
import tarfile
import tempfile
from pathlib import Path
import pandas as pd
import anndata as ad
from scipy.io import mmread

# Resolved from this file's location so paths work regardless of cwd (scripts, notebooks/)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "output"

"""
See AnnData documentation at  
  https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html#conversion-to-dataframes
and
  https://www.sc-best-practices.org/introduction/fundamental-data-structures-and-frameworks/
"""

###
# Reading BROAD Institute scRNA files (matrix.mtx, barcodes.tsv, features.tsv).
# Files may carry a sample prefix (e.g. Sample_3444_matrix.mtx) and may be gzipped.
###
def _find_file(data_dir, suffixes):
    """Return the first file in data_dir whose name ends with one of `suffixes`."""
    for name in sorted(os.listdir(data_dir)):
        if name.endswith(tuple(suffixes)):
            return os.path.join(data_dir, name)
    raise FileNotFoundError(f"None of {suffixes} found in {data_dir}")


def load_broad_dataset(data_dir):
    """
    Reads an already-extracted 10x-style trio (matrix / barcodes / features|genes)
    from `data_dir` into an AnnData (cells x genes, float32 CSR).
    """
    matrix_path = _find_file(data_dir, ["matrix.mtx", "matrix.mtx.gz"])
    barcodes_path = _find_file(data_dir, ["barcodes.tsv", "barcodes.tsv.gz"])
    genes_path = _find_file(data_dir, ["features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"])

    # mtx is genes x cells; AnnData wants cells x genes. Float for downstream normalization.
    X = mmread(matrix_path).tocsr().T.astype("float32")

    barcodes = pd.read_csv(barcodes_path, header=None, sep='\t', names=['barcode'])
    genes = pd.read_csv(genes_path, header=None, sep='\t')

    # features.tsv has 3 cols (id, symbol, type); old genes.tsv has 2; some have only ids
    genes.columns = ['gene_ids', 'gene_symbols', 'feature_type'][:genes.shape[1]]

    barcodes.index = barcodes['barcode'].astype(str)
    if 'gene_symbols' in genes.columns:
        genes.index = ad.utils.make_index_unique(pd.Index(genes['gene_symbols'].astype(str)))
    else:
        genes.index = genes['gene_ids'].astype(str)

    if X.shape != (len(barcodes), len(genes)):
        raise ValueError(f"Shape mismatch: matrix {X.shape} vs {len(barcodes)} barcodes x {len(genes)} genes")

    return ad.AnnData(X=X, obs=barcodes, var=genes)


def _sample_name(sample_id):
    """3444 / '3444' / 'Sample_3444' -> 'Sample_3444'."""
    s = str(sample_id)
    return s if s.startswith("Sample_") else f"Sample_{s}"


def load_broad_metadata(metadata_path=DATA_DIR / "metadata.json"):
    """Returns metadata.json as a DataFrame indexed by sample name ('Sample_3444')."""
    meta = pd.read_json(metadata_path, dtype=False)
    meta["age"] = pd.to_numeric(meta["age"])
    meta.index = meta["sample_id"].map(_sample_name)
    return meta


def load_broad_sample(sample_id, data_dir=DATA_DIR, metadata_path=None):
    """
    Unpacks data_dir/Sample_<id>.tar(.gz) into a temp dir, reads it, and returns an AnnData.

    Sample metadata (age, sex, tissue, disease_status, ...) is added to every cell in
    `adata.obs` and the full row is kept in `adata.uns['metadata']`.
    The temp dir is deleted afterwards; only the AnnData is kept.
    """
    name = _sample_name(sample_id)
    archive = next((os.path.join(data_dir, name + ext) for ext in (".tar.gz", ".tar")
                    if os.path.exists(os.path.join(data_dir, name + ext))), None)
    if archive is None:
        raise FileNotFoundError(f"No {name}.tar.gz / {name}.tar in {data_dir}")

    with tempfile.TemporaryDirectory() as tmp:
        with tarfile.open(archive) as tar:
            tar.extractall(tmp, filter="data")  # 'data' filter blocks unsafe paths
        # archives contain a single Sample_X/ folder; fall back to tmp if flat
        sub = os.path.join(tmp, name)
        adata = load_broad_dataset(sub if os.path.isdir(sub) else tmp)

    metadata_path = metadata_path or os.path.join(data_dir, "metadata.json")
    meta = load_broad_metadata(metadata_path)
    if name not in meta.index:
        raise KeyError(f"{name} not found in {metadata_path}")
    row = meta.loc[name]

    adata.obs["sample"] = name
    for col in ["age", "sex", "tissue", "disease_status", "study_id", "format"]:
        adata.obs[col] = row[col]
    adata.uns["metadata"] = row.to_dict()
    return adata


def load_broad_samples(sample_ids=None, data_dir=DATA_DIR, metadata_path=None):
    """Loads several samples (default: every sample in metadata.json) -> {'Sample_3444': AnnData}."""
    if sample_ids is None:
        meta = load_broad_metadata(metadata_path or os.path.join(data_dir, "metadata.json"))
        sample_ids = list(meta.index)
    return {_sample_name(i): load_broad_sample(i, data_dir, metadata_path) for i in sample_ids}


###
# Reading .h5ad scRNA type files.
###
def load_h5ad_dataset(data_dir):
    # Read the file into an AnnData object
    adata = ad.read_h5ad(data_dir)
    return adata




if __name__ == "__main__":
    adata = load_broad_sample(3444)
    print(adata)
    print(adata.obs.head())
