#NAME: convert_to_h5ad.py
# ST
"""
Converts raw data/Sample_<id>.tar(.gz) archives to data/Sample_<id>.h5ad.

Usage (from the project root):
  python src/01_preprocess/convert_to_h5ad.py                    # every sample in metadata.json
  python src/01_preprocess/convert_to_h5ad.py --samples 3444 1931
  python src/01_preprocess/convert_to_h5ad.py --overwrite
"""

# IMPORTS
import argparse
import anndata as ad
from cs223_scRNAdataset_read_funs import DATA_DIR, _sample_name, load_broad_metadata, load_broad_sample


def h5ad_path(sample_id):
    return DATA_DIR / f"{_sample_name(sample_id)}.h5ad"


def convert(sample_ids=None, overwrite=False):
    """Writes one .h5ad per sample; existing files are skipped unless overwrite=True."""
    if sample_ids is None:
        sample_ids = list(load_broad_metadata().index)

    for sample_id in sample_ids:
        out = h5ad_path(sample_id)
        if out.exists() and not overwrite:
            print(f"skip   {out.name} (exists)")
            continue
        adata = load_broad_sample(sample_id)
        adata.write_h5ad(out, compression="gzip")
        print(f"wrote  {out.name}: {adata.n_obs} cells x {adata.n_vars} genes")


def load_h5ad_sample(sample_id):
    """Reads a converted sample back: 3444 / 'Sample_3444' -> AnnData."""
    return ad.read_h5ad(h5ad_path(sample_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", nargs="+", help="sample ids (default: all in metadata.json)")
    parser.add_argument("--overwrite", action="store_true", help="rewrite existing .h5ad files")
    args = parser.parse_args()
    convert(args.samples, args.overwrite)
