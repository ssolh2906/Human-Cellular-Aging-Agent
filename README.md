# Human-Cellular-Aging-Agent

## Setup: Google Drive credentials

To download data from a restricted Google Drive, you need OAuth credentials.

1. Go to [Google Cloud Console](https://console.cloud.google.com) and create a new project
2. **APIs & Services → Library** → search `Google Drive API` → Enable
3. **APIs & Services → Credentials** → Create `OAuth 2.0 Client ID` (type: **Desktop app**)
4. Download `credentials.json` → save it as `config/credentials.json`

> `config/credentials.json` and `config/token.json` are listed in `.gitignore` and will not be committed.

```
Human-Cellular-Aging-Agent/
├── data/                  # raw h5ad
├── output/                # processed h5ad
├── analysis/              # scBubbleTree
├── analysis_results/      # project analysis
├── src/
│   ├── 01_preprocess/
│   ├── 02_bubbletree/
│   ├── 03_grn/
│   ├── 04_pseudotime/
│   ├── 05_perturbation/
│   └── 06_validation/
├── notebooks/
├── docs/
├── logs/
└── requirements.txt
```
