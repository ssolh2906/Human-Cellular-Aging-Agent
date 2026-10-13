# Human Cellular Aging Agent — Task Overview

## Background

Cells age by losing epigenetic marks, degrading function and regenerative capacity. The goal of this project is to find **transcription factors (TFs)** that can partially reprogram old cells back to a youthful state — without erasing cell identity or causing tumor risk.

Key concepts:
- **Partial epigenetic reprogramming**: reset the epigenome using TF combinations
- **AI discovery engine**: use in-silico experiments to narrow the massive TF search space
- **LNP-mRNA delivery**: temporarily express rejuvenating TFs via lipid nanoparticles

---

## Objective

Build a computational pipeline using **scRNA-seq data** to:
1. Characterize aging at the single-cell level
2. Infer which transcription factors regulate aging genes
3. Simulate which TF perturbations reverse the aging signature

---

## Pipeline Steps

### Step 1 — Define the Aging Transcriptomic Signature

**Goal:** Establish what "young" vs "old" looks like at single-cell resolution.

**Data setup:**
- Download datasets from the "Cell Aging" Google Drive folder
- Each `Sample_<#>.tar` unpacks to 3 files: `_matrix.mtx`, `_feature.tsv`, `_barcodes.tsv`
- `metadata.json` maps each sample to young/elderly age group
- Convert all `.mtx`/`.tsv` trios to `.h5ad` format (sample code provided)
- Save `.h5ad` files to `./data/`

**Analysis (via Claude Sonnet):**
- Read each `.h5ad` and run standard scRNA-seq pipeline: QC → normalization → clustering
- Save processed outputs to `./output/` with a modified filename

**scBubbleTree:**
- Build a scBubbleTree for each processed `.h5ad` output
- Use trees to:
  - **Differential Expression (DE):** compare young vs old within the same cell type
  - **Age-associated markers:** identify up-regulated genes (senescence, inflammation, dysfunction) and down-regulated genes (metabolism, proteostasis)

---

### Step 2 — Infer Gene Regulatory Networks (GRNs)

**Goal:** Map which TFs control age-associated genes.

- **SCENIC**: infer active TFs and their target gene sets (regulons)
- **Network topology**: identify "hub" TFs that regulate large modules of aging genes

---

### Step 3 — Model Cellular Trajectories & Pseudotime

**Goal:** Treat aging as a continuous process, not a binary state.

- **Pseudotime analysis** (Monocle3 or Palantir): order cells along a young → old trajectory
- **Reprogramming trajectories**: if Yamanaka factor data is available, map which proteins drive cells backward along the aging axis

---

### Step 4 — In Silico Perturbation Screening

**Goal:** Computationally predict which TF perturbations reverse aging.

- **CellOracle**: simulate TF over-expression or knockout; predict network-wide gene expression shifts
- **Foundation models**: use scGPT or GEARS to predict transcriptomic outcomes of gene perturbations
- **Reversal scoring**: rank perturbations by how closely the predicted old-cell profile matches an actual young-cell profile

---

### Step 5 — Prioritize and Validate Candidates

**Goal:** Filter computational hits to a shortlist for wet-lab follow-up.

- **Druggability check**: cross-reference candidates with DGIdb to find targetable proteins
- **Experimental validation plan**: CRISPRa, RNAi, or small molecules in aged cell cultures; measure rejuvenation using epigenetic or transcriptomic aging clocks

---

## Key Tools

| Task | Tool |
|---|---|
| scRNA-seq analysis | Scanpy / Claude Sonnet |
| Visualization | scBubbleTree |
| GRN inference | SCENIC |
| Pseudotime | Monocle3 / Palantir |
| Network perturbation | CellOracle |
| DL perturbation prediction | scGPT / GEARS |
| Druggability lookup | DGIdb |

---

## File Conventions

| Path | Contents |
|---|---|
| `data/` | Input `.h5ad` files (converted from raw `.mtx`/`.tsv`) |
| `output/` | Processed `.h5ad` files after QC + normalization |
| `results/` | Figures, reports, candidate lists |
| `src/01_preprocess/` | Data loading, QC, normalization scripts |
| `src/02_bubbletree/` | scBubbleTree construction + DE analysis |
| `src/03_grn/` | SCENIC GRN inference |
| `src/04_pseudotime/` | Monocle3 / Palantir trajectory analysis |
| `src/05_perturbation/` | CellOracle, scGPT, GEARS perturbation scripts |
| `src/06_validation/` | Candidate prioritization and DGIdb lookup |
| `notebooks/` | Exploratory Jupyter notebooks |
