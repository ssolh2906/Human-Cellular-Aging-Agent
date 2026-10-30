#ST
"""
Ensembl ID -> gene symbol / chromosome table. The raw data only has ENSG IDs, so MT- genes
and marker genes can't be found by name without it.

The table (data/gene_map.tsv, gitignored) is built once from the Ensembl release-98 GTF
(the release the 10x GRCh38-2020 reference is based on). If it is missing, run:
  python src/01_preprocess/gene_map.py          # asks y/N, then downloads (~47 MB, keeps ~2 MB)

Use in code:
  from gene_map import add_gene_info
  add_gene_info(adata)    # adds var columns: symbol, chrom, is_mt
"""

# IMPORTS
import gzip
import re
import sys
import urllib.request

import pandas as pd

from cs223_scRNAdataset_read_funs import DATA_DIR

GENE_MAP_PATH = DATA_DIR / "gene_map.tsv"
GTF_URL = "https://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.chr.gtf.gz"
DOWNLOAD_NOTE = f"{GTF_URL} (~47 MB; only the gene rows, ~2 MB, are kept in {GENE_MAP_PATH.name})"


def _attr(attrs, key):
    m = re.search(rf'{key} "([^"]+)"', attrs)
    return m.group(1) if m else ""


def download_gene_map():
    """Streams the GTF and keeps only `gene` rows; the big file is never stored."""
    rows = []
    with urllib.request.urlopen(GTF_URL, timeout=120) as resp, gzip.open(resp, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            chrom, _, feature, *_, attrs = line.rstrip("\n").split("\t")
            if feature == "gene":
                rows.append((_attr(attrs, "gene_id"), _attr(attrs, "gene_name"), chrom, _attr(attrs, "gene_biotype")))
    table = pd.DataFrame(rows, columns=["gene_id", "symbol", "chrom", "biotype"])
    DATA_DIR.mkdir(exist_ok=True)
    tmp = GENE_MAP_PATH.with_suffix(".tmp")   # write then rename, so a failed download leaves no half file
    table.to_csv(tmp, sep="\t", index=False)
    tmp.replace(GENE_MAP_PATH)
    return table


def load_gene_map():
    if not GENE_MAP_PATH.exists():
        raise FileNotFoundError(
            f"{GENE_MAP_PATH} is missing. Download it with the fetch_gene_map tool "
            "or `python src/01_preprocess/gene_map.py`."
        )
    return pd.read_csv(GENE_MAP_PATH, sep="\t", keep_default_na=False).set_index("gene_id")


def add_gene_info(adata):
    """Adds var['symbol'], var['chrom'], var['is_mt'] by Ensembl ID (var_names). Unmatched IDs keep their ID as symbol."""
    table = load_gene_map()
    ids = adata.var_names.str.replace(r"\.\d+$", "", regex=True)   # drop version suffix if present
    matched = ids.isin(table.index)
    info = table.reindex(ids)
    adata.var["symbol"] = [s if ok and s else i for s, ok, i in zip(info["symbol"], matched, ids)]
    adata.var["chrom"] = info["chrom"].fillna("").to_numpy()
    adata.var["is_mt"] = (adata.var["chrom"] == "MT").to_numpy()
    print(f"gene map: {matched.sum()}/{len(ids)} genes matched ({matched.mean():.1%}), {adata.var['is_mt'].sum()} MT genes")
    return adata


if __name__ == "__main__":
    if GENE_MAP_PATH.exists():
        sys.exit(f"{GENE_MAP_PATH} already exists; delete it to rebuild.")
    print(f"Will download: {DOWNLOAD_NOTE}")
    if input("Download? [y/N] ").strip().lower() != "y":
        sys.exit("cancelled")
    table = download_gene_map()
    print(f"wrote {GENE_MAP_PATH.name}: {len(table)} genes")
