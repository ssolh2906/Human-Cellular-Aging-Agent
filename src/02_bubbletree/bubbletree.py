#ST
"""
scBubbleTree in pure Python (not the R package), rebuilt from the course examples.

A tree is built from the clusters (e.g. Leiden) of ONE processed AnnData:
  leaf bubbles = clusters, internal nodes = ward merges of the cluster centroids.
Every node carries a centroid in PCA space (for distances / A*) and a mean expression
vector over all genes (for cross-tree comparison: PCA spaces of different samples are not comparable).

  tree = build_bubbletree(adata, name="Sample_3444")
  save_tree(tree)                              # -> analysis/bubbletree_<date>_<time>/Sample_3444_bubbletree.*
  trees = load_trees(latest_run_dir())         # later: all samples of a run, for the joint analysis
  astar_path(tree, "0", "3")
  find_profile(trees, {"CDKN1A": "high", "MKI67": "low"})
"""

# IMPORTS
import heapq
import json
import re
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import cdist

# Resolved from this file's location so paths work regardless of cwd (scripts, notebooks/, agent kernel)
PROJECT_ROOT = Path(__file__).resolve().parents[2]
ANALYSIS_DIR = PROJECT_ROOT / "analysis"
TREE_SUFFIX = "_bubbletree.h5ad"
DEFAULT_NUMERIC = ("total_counts", "n_genes_by_counts", "pct_counts_mt")
# node-table columns that are structure, not per-bubble statistics
STRUCT_COLS = {"is_leaf", "height", "parent", "left", "right"}


# ---------------------------------------------------------------- tree object
class BubbleTree:
    """
    nodes:  one row per node. Rows 0..k-1 are the leaf bubbles, k..2k-2 the internal nodes
            (same numbering as scipy's linkage matrix Z, so Z row i is node k+i).
    coords: node centroids in PCA space (cell-count weighted for internal nodes)
    expr:   node x gene mean expression (same weighting), genes = var.index
    """

    def __init__(self, name, nodes, coords, expr, var, Z, params):
        self.name, self.nodes, self.coords, self.expr = name, nodes, coords, expr
        self.var, self.Z, self.params = var, Z, params
        self.n_leaves = int(nodes["is_leaf"].sum())

    @property
    def leaves(self):
        return self.nodes.iloc[: self.n_leaves]

    @property
    def stat_cols(self):
        return [c for c in self.nodes.columns if c not in STRUCT_COLS]

    def node_index(self, node):
        node = str(node)
        if node not in self.nodes.index:
            raise ValueError(f"unknown node {node!r}; leaf bubbles are {list(self.leaves.index)}")
        return self.nodes.index.get_loc(node)

    def find_gene(self, gene):
        """Column of `expr` for an Ensembl ID or (case-insensitive) gene symbol; None if absent."""
        if gene in self.var.index:
            return self.var.index.get_loc(gene)
        if "symbol" in self.var:
            hits = np.flatnonzero(self.var["symbol"].astype(str).str.upper().to_numpy() == gene.upper())
            if len(hits):
                return int(hits[0])
        return None

    def gene_values(self, gene):
        j = self.find_gene(gene)
        if j is None:
            raise ValueError(f"gene {gene!r} not found in tree {self.name} (Ensembl ID or symbol)")
        return self.expr[:, j]


# ---------------------------------------------------------------- build
def _natural(label):
    """'2' before '10', numbers before names."""
    return (0, int(label), "") if label.isdigit() else (1, 0, label)


def _warn_if_not_lognorm(X):
    x = X[:1000]
    vals = x.data if sp.issparse(x) else np.asarray(x).ravel()
    if vals.size and vals.min() < 0:
        warnings.warn("X has negative values (scaled data?); bubble expression should be log-normalized. "
                      "Pass layer=... with log-normalized values.")
    elif vals.size and np.all(vals == np.round(vals)):
        warnings.warn("X looks like raw counts; bubble expression should be log-normalized. "
                      "Pass layer=... with log-normalized values.")


def _propagate(leaf_vals, sizes, Z):
    """Leaf values (k x d) -> all nodes (2k-1 x d): each merge is the cell-count weighted mean of its children."""
    k = len(sizes)
    vals = np.zeros((2 * k - 1, leaf_vals.shape[1]), dtype=leaf_vals.dtype)
    vals[:k] = leaf_vals
    n = np.concatenate([sizes, np.zeros(k - 1)]).astype(float)
    for i, (a, b) in enumerate(Z[:, :2].astype(int)):
        n[k + i] = n[a] + n[b]
        vals[k + i] = (n[a] * vals[a] + n[b] * vals[b]) / n[k + i]
    return vals, n


def build_bubbletree(adata, name, cluster_key="leiden", embed_key="X_pca", n_pcs=None, layer=None,
                     group_cols=(), numeric_cols=DEFAULT_NUMERIC):
    """
    adata: processed, with clusters in obs[cluster_key], a PCA in obsm[embed_key] and log-normalized X
           (or the log-normalized matrix in `layer`).
    group_cols: obs columns whose per-bubble composition is stored as frac__<col>__<level> (e.g. age_group).
    numeric_cols: obs columns whose per-bubble mean is stored as mean_<col> (skipped if absent).
    """
    if cluster_key not in adata.obs:
        raise ValueError(f"obs has no {cluster_key!r}")
    if embed_key not in adata.obsm:
        raise ValueError(f"obsm has no {embed_key!r}")

    raw_labels = adata.obs[cluster_key].astype(str).to_numpy()
    labels = sorted(set(raw_labels), key=_natural)
    k, n = len(labels), adata.n_obs
    if k < 2:
        raise ValueError("need at least 2 clusters to build a tree")
    code = pd.Categorical(raw_labels, categories=labels).codes.astype(int)
    sizes = np.bincount(code, minlength=k)

    # Row-normalised cluster indicator: M @ A is the per-cluster mean of A (works for sparse A)
    M = sp.csr_matrix((1.0 / sizes[code], (code, np.arange(n))), shape=(k, n))
    X = adata.layers[layer] if layer else adata.X
    _warn_if_not_lognorm(X)
    centroids = np.asarray(M @ np.asarray(adata.obsm[embed_key])[:, :n_pcs])
    leaf_expr = M @ X
    leaf_expr = np.asarray(leaf_expr.todense() if sp.issparse(leaf_expr) else leaf_expr, dtype=np.float32)

    stats = pd.DataFrame(index=labels)
    for col in numeric_cols:
        if col in adata.obs:
            stats[f"mean_{col}"] = pd.Series(adata.obs[col].to_numpy()).groupby(code).mean().to_numpy()
    for col in group_cols:
        comp = pd.crosstab(code, adata.obs[col].astype(str).to_numpy(), normalize="index")
        for level in comp.columns:
            stats[f"frac__{col}__{level}"] = comp[level].to_numpy()

    # Ward on the centroids, as in the course example (scBubbletree also clusters cluster-level profiles)
    Z = linkage(centroids, method="ward")
    coords, n_all = _propagate(centroids, sizes, Z)
    expr, _ = _propagate(leaf_expr, sizes, Z)
    stat_vals, _ = _propagate(stats.to_numpy(dtype=float), sizes, Z)

    names = labels + [f"node_{i}" for i in range(k - 1)]
    parent = [""] * (2 * k - 1)
    left, right = [""] * (2 * k - 1), [""] * (2 * k - 1)
    for i, (a, b) in enumerate(Z[:, :2].astype(int)):
        parent[a] = parent[b] = names[k + i]
        left[k + i], right[k + i] = names[a], names[b]
    nodes = pd.DataFrame({
        "is_leaf": [True] * k + [False] * (k - 1),
        "n_cells": n_all.astype(int),
        "frac": n_all / n,
        "height": np.concatenate([np.zeros(k), Z[:, 2]]),   # ward merge distance
        "parent": parent, "left": left, "right": right,
    }, index=names)
    for j, col in enumerate(stats.columns):
        nodes[col] = stat_vals[:, j]

    var = pd.DataFrame(index=adata.var_names.astype(str))
    if "symbol" in adata.var:
        var["symbol"] = adata.var["symbol"].astype(str).to_numpy()   # lets find_gene accept symbols
    params = {"cluster_key": cluster_key, "embed_key": embed_key, "n_pcs": n_pcs, "layer": layer,
              "n_cells": n, "n_clusters": k}
    return BubbleTree(name, nodes, coords, expr, var, Z, params)


# ---------------------------------------------------------------- plot
def plot_bubbletree(tree, color_by=None, highlight=(), save_path=None, cmap="viridis"):
    """
    Dendrogram with one bubble per leaf: area = cell count, color = `color_by`
    (a node-table column such as mean_n_genes_by_counts / frac__age_group__old, or a gene).
    highlight: leaf names drawn with a red outline (e.g. an A* path). Saves a PNG if save_path is given.
    """
    k = tree.n_leaves
    if color_by is None:
        color_by = next((c for c in tree.stat_cols if c.startswith("mean_")), "n_cells")
    values = (tree.nodes[color_by].to_numpy() if color_by in tree.nodes else tree.gene_values(color_by))[:k]

    fig, ax = plt.subplots(figsize=(max(7, 1.1 * k), 6))
    dend = dendrogram(tree.Z, labels=list(tree.leaves.index), ax=ax, link_color_func=lambda _: "#4A5568")
    order = dend["leaves"]
    xs = 5 + 10 * np.arange(k)   # scipy places the i-th leaf at x = 5 + 10*i
    leaves = tree.leaves.iloc[order]

    n_cells = leaves["n_cells"].to_numpy()
    span = (n_cells.max() - n_cells.min()) or 1
    sizes = 400 + (n_cells - n_cells.min()) / span * 3100
    edge = ["red" if name in set(map(str, highlight)) else "black" for name in leaves.index]
    sc = ax.scatter(xs, np.zeros(k), s=sizes, c=values[order], cmap=cmap, edgecolors=edge,
                    linewidths=[3 if e == "red" else 1.5 for e in edge], zorder=3, clip_on=False)
    for x, name, n_c, s in zip(xs, leaves.index, n_cells, sizes):
        # offset in points (bubble radius + gap), so the label clears the bubble whatever the tree height is
        ax.annotate(f"C{name}\nn={n_c}", (x, 0), xytext=(0, -(np.sqrt(s / np.pi) + 6)), textcoords="offset points",
                    ha="center", va="top", fontsize=9, fontweight="bold")

    fig.colorbar(sc, ax=ax, orientation="horizontal", pad=0.18, shrink=0.6).set_label(color_by)
    ax.set_title(f"scBubbleTree: {tree.name}", fontsize=14, pad=15)
    ax.set_ylabel("Ward distance between bubbles", fontsize=12)
    ax.set_xticks([])
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    fig.tight_layout()

    if save_path is None:
        return fig
    fig.savefig(save_path, dpi=120)
    plt.close(fig)
    return Path(save_path)


# ---------------------------------------------------------------- A*
def _build_graph(tree, knn):
    """
    Tree edges (parent-child) plus kNN bridge edges between leaves. In a pure tree there is only one path,
    so A* would have nothing to decide. Edge weight = Euclidean distance between the node centroids.
    """
    k = tree.n_leaves
    adj = {i: {} for i in range(2 * k - 1)}

    def link(a, b):
        adj[a][b] = adj[b][a] = float(np.linalg.norm(tree.coords[a] - tree.coords[b]))

    for i, (a, b) in enumerate(tree.Z[:, :2].astype(int)):
        link(k + i, a)
        link(k + i, b)
    d = cdist(tree.coords[:k], tree.coords[:k])
    np.fill_diagonal(d, np.inf)
    for i in range(k):
        for j in np.argsort(d[i])[: min(knn, k - 1)]:
            link(i, int(j))
    return adj


def astar_path(tree, start, goal, knn=2):
    """
    A* between two nodes of one tree, g = path length, h = straight Euclidean distance to the goal in PCA space.
    The course example used edge weights that could be shorter than the straight line, which breaks
    admissibility; with Euclidean edge weights h is admissible and consistent, so the path is optimal.
    """
    s, g = tree.node_index(start), tree.node_index(goal)
    adj = _build_graph(tree, knn)
    h = lambda i: float(np.linalg.norm(tree.coords[i] - tree.coords[g]))

    heap, counter = [(h(s), 0, s)], 0   # counter breaks f-score ties without comparing nodes
    g_score, came_from, closed = {s: 0.0}, {}, set()
    while heap:
        _, _, cur = heapq.heappop(heap)
        if cur in closed:
            continue   # stale entry: a cheaper route to this node was pushed later
        if cur == g:
            path = [cur]
            while path[-1] in came_from:
                path.append(came_from[path[-1]])
            return {"path": [tree.nodes.index[i] for i in path[::-1]], "cost": g_score[g],
                    "n_expanded": len(closed)}
        closed.add(cur)
        for nb, w in adj[cur].items():
            tentative = g_score[cur] + w
            if tentative < g_score.get(nb, np.inf):
                g_score[nb], came_from[nb] = tentative, cur
                counter += 1
                heapq.heappush(heap, (tentative + h(nb), counter, nb))
    raise RuntimeError("no path found")   # unreachable: tree edges keep the graph connected


# ---------------------------------------------------------------- cross-tree search
def find_profile(trees, target, top=10, min_cells=10):
    """
    Which bubbles, over all given trees, best match a gene profile? target = {gene: "high" | "low"}.
    Each gene is z-scored over the POOLED leaf bubbles of all trees, so scores are comparable across trees;
    score = mean over genes of (+z for high, -z for low). Higher = closer to the profile.
    Replaces the course example's cross-tree A*: with tens of bubbles per tree an exact scan finds the best
    match directly, and that example's heuristic was not admissible.
    min_cells drops tiny bubbles whose mean expression is noise.
    """
    trees = list(trees.values()) if isinstance(trees, dict) else [trees] if isinstance(trees, BubbleTree) else list(trees)
    if any(v not in ("high", "low") for v in target.values()):
        raise ValueError('target values must be "high" or "low"')
    genes = list(target)

    parts, mats = [], []
    for t in trees:
        k = t.n_leaves
        keep = (t.leaves["n_cells"] >= min_cells).to_numpy()
        E = np.full((k, len(genes)), np.nan)
        for j, gene in enumerate(genes):
            col = t.find_gene(gene)
            if col is not None:
                E[:, j] = t.expr[:k, col]
        part = t.leaves[["n_cells", "frac"] + [c for c in t.stat_cols if c not in ("n_cells", "frac")]]
        part = part.copy()
        part.insert(0, "bubble", part.index)
        part.insert(0, "tree", t.name)
        parts.append(part[keep])
        mats.append(E[keep])
    table, E = pd.concat(parts, ignore_index=True), np.vstack(mats)

    absent = [g for j, g in enumerate(genes) if np.isnan(E[:, j]).all()]
    if absent:
        raise ValueError(f"genes not found in any tree: {absent}")
    if np.isnan(E).any():
        warnings.warn("some genes are missing in some trees; those trees are scored on the remaining genes")

    mu, sd = np.nanmean(E, axis=0), np.nanstd(E, axis=0)
    z = (E - mu) / np.where(sd == 0, 1, sd)
    sign = np.array([1.0 if target[g] == "high" else -1.0 for g in genes])
    table.insert(2, "score", np.nanmean(z * sign, axis=1))
    for j, gene in enumerate(genes):
        table[f"expr_{gene}"] = E[:, j]
    return table.sort_values("score", ascending=False).head(top).reset_index(drop=True)


# ---------------------------------------------------------------- save / load
_RUN_DIR = None


def make_run_dir(base=None):
    """New analysis/bubbletree_<YYYYMMDD_HHMMSS>/ directory."""
    d = Path(base or ANALYSIS_DIR) / f"bubbletree_{datetime.now():%Y%m%d_%H%M%S}"
    d.mkdir(parents=True, exist_ok=True)
    return d


def get_run_dir():
    """One dated directory per Python process, so every tree saved in one agent run lands together."""
    global _RUN_DIR
    if _RUN_DIR is None:
        _RUN_DIR = make_run_dir()
    return _RUN_DIR


def latest_run_dir(base=None):
    runs = sorted(Path(base or ANALYSIS_DIR).glob("bubbletree_*"))
    if not runs:
        raise FileNotFoundError("no bubbletree_* run directory yet")
    return runs[-1]


def save_tree(tree, run_dir=None, plot=True):
    """
    Writes <name>_bubbletree.h5ad (everything needed to reload the tree), <name>_bubbles.csv (leaf table to read)
    and <name>_bubbletree.png into run_dir (default: this process's dated directory). Returns the paths.
    """
    run_dir = Path(run_dir) if run_dir else get_run_dir()
    run_dir.mkdir(parents=True, exist_ok=True)
    stem = re.sub(r"[^\w.-]", "_", tree.name)

    # nodes x genes AnnData: stays an ordinary h5ad, so trees of different samples can be concatenated by gene
    store = ad.AnnData(X=tree.expr.astype(np.float32), obs=tree.nodes.copy(), var=tree.var.copy())
    store.obsm["coords"] = tree.coords
    store.uns["linkage"] = tree.Z
    store.uns["name"] = tree.name
    store.uns["params_json"] = json.dumps(tree.params)
    h5ad = run_dir / f"{stem}{TREE_SUFFIX}"
    store.write_h5ad(h5ad)

    csv = run_dir / f"{stem}_bubbles.csv"
    tree.leaves.drop(columns=["is_leaf", "parent", "left", "right"]).to_csv(csv, index_label="bubble")
    out = {"h5ad": h5ad, "csv": csv}
    if plot:
        out["png"] = plot_bubbletree(tree, save_path=run_dir / f"{stem}_bubbletree.png")
    return out


def load_tree(path):
    a = ad.read_h5ad(path)
    nodes = a.obs.copy()
    for col in ("parent", "left", "right"):
        nodes[col] = nodes[col].astype(str).replace("nan", "")
    nodes["is_leaf"] = nodes["is_leaf"].astype(bool)
    expr = a.X.toarray() if sp.issparse(a.X) else np.asarray(a.X)
    return BubbleTree(str(a.uns["name"]), nodes, np.asarray(a.obsm["coords"]), expr, a.var.copy(),
                      np.asarray(a.uns["linkage"]), json.loads(a.uns["params_json"]))


def load_trees(run_dir):
    """All trees of one run directory as {name: BubbleTree}."""
    trees = (load_tree(p) for p in sorted(Path(run_dir).glob(f"*{TREE_SUFFIX}")))
    return {t.name: t for t in trees}
