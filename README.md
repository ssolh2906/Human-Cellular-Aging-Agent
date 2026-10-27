# Human-Cellular-Aging-Agent

## Setup: Google Drive credentials

To download data from a restricted Google Drive, you need OAuth credentials.

1. Go to [Google Cloud Console](https://console.cloud.google.com) and create a new project
2. **APIs & Services → Library** → search `Google Drive API` → Enable
3. **APIs & Services → Credentials** → Create `OAuth 2.0 Client ID` (type: **Desktop app**)
4. Download `credentials.json` → save it as `config/credentials.json`

> `config/credentials.json` and `config/token.json` are listed in `.gitignore` and will not be committed.

## Setup: Claude API key

The analysis agent (`src/agent/`) calls the Claude API.

1. Get an API key from the [Claude Console](https://platform.claude.com) (**API Keys → Create Key**)
2. Create a `.env` file in the project root:
   ```
   ANTHROPIC_API_KEY=sk-ant-...
   ```
3. The agent loads it with `python-dotenv`

> `.env` is listed in `.gitignore` and will not be committed. Don't put the key in a IDE run configuration's environment variables: a shared run configuration (`.run/`) stores it in plain text.

<details>
<summary><b>PyCharm run configuration</b></summary>

Run → Edit Configurations… → *Edit configuration templates…* → Python, then set
- **Python interpreter:** the project `.venv`
- **Working directory:** the project root
- **Paths to ".env" files:** `.env`

Every new Python run configuration inherits these.

</details>

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
