# Agent Coding Conventions

## Pattern: ReAct (Reason → Act → Observe)

This project uses a **tool-use agentic loop** via the Anthropic API.
Claude reasons about the next step, calls a tool (e.g. `execute_python`), observes the result, and repeats until the task is complete.

```
Thought  →  decide what analysis to run
Act      →  execute_python(code)
Observe  →  read stdout / errors / file output
Thought  →  adjust and continue, or declare done
```

---

## Project Layout

```
src/<step>/
├── agent.py        ← agent loop (hand-written)
├── tools.py        ← tool definitions (hand-written)
└── pipeline/
    ├── Sample_1.py ← final confirmed script (Claude-generated, saved once)
    └── Sample_2.py

logs/               ← full run logs (gitignored)
└── Sample_1_run_20261013_142305.json
```

---

## What to Save

### Save: final pipeline script only

Save the generated script **only when Claude declares the analysis complete and correct** — not on every attempt.

- Location: `src/<step>/pipeline/<sample_name>.py`
- One file per sample per pipeline step
- The file must be runnable standalone (no agent loop required)

```python
# In tools.py — save only on success signal
def save_final_script(code: str, sample: str, step: str):
    path = f"src/{step}/pipeline/{sample}.py"
    with open(path, "w") as f:
        f.write(code)
```

### Save: full conversation log on every run

Save the entire `messages` array as JSON after each run, regardless of success or failure.

- Location: `logs/<sample>_run_<timestamp>.json`
- Contains every attempt, error, and reasoning step
- Used for debugging — never committed to git

```python
import json, datetime

def save_log(messages: list, sample: str):
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = f"logs/{sample}_run_{ts}.json"
    with open(path, "w") as f:
        json.dump(messages, f, indent=2)
```

---

## What NOT to Save

| What | Why |
|---|---|
| Intermediate code attempts | Noise — the log already has them |
| Partial / erroring scripts | Misleading as standalone files |
| Raw API response objects | Redundant with the messages log |

---

## File Naming

| Type | Pattern | Example |
|---|---|---|
| Final pipeline script | `<Sample_name>.py` | `Sample_1.py` |
| Run log | `<Sample_name>_run_<YYYYMMDD_HHMMSS>.json` | `Sample_1_run_20261013_142305.json` |
| Processed data output | `<Sample_name>_processed.h5ad` | `Sample_1_processed.h5ad` |

---

## Git Policy

| Path | Tracked |
|---|---|
| `src/<step>/agent.py` | Yes |
| `src/<step>/tools.py` | Yes |
| `src/<step>/pipeline/*.py` | Yes — final scripts only |
| `logs/` | **No** — gitignored |
| `output/*.h5ad` | **No** — gitignored (large files) |
