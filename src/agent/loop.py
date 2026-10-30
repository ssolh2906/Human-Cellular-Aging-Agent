#ST
"""
Hand-written ReAct loop (reason -> tool call -> observe) on the Anthropic Messages API.
Manual loop, not the Tool Runner: we must keep the full `messages` log of every run.
"""

# IMPORTS
import datetime
import json
from pathlib import Path

import anthropic

from kernel import PythonKernel
from tools import TOOLS, run_tool

PROJECT_ROOT = Path(__file__).resolve().parents[2]
LOG_DIR = PROJECT_ROOT / "logs"
PROMPT_DIR = Path(__file__).resolve().parent / "prompts"

MODEL = "claude-sonnet-5-5"
MAX_TOKENS = 32000
# Sonnet 5.5 list prices, USD per 1M tokens (cache write assumed 1.25x input) - estimate only
PRICE = {"input": 2.0, "output": 10.0, "cache_read": 0.2, "cache_write": 2.5}


def load_prompt(name="system"):
    """Reads prompts/<name>.md, dropping the '#ST' marker line."""
    lines = (PROMPT_DIR / f"{name}.md").read_text().splitlines()
    return "\n".join(l for l in lines if l.strip() != "#ST").strip()


def _drop_image_data(obj):
    """Copy of the log with base64 image bytes removed (they would bloat every log by MBs)."""
    if isinstance(obj, dict):
        if obj.get("type") == "base64" and "data" in obj:
            return {**obj, "data": f"<{len(obj['data'])} base64 chars omitted>"}
        return {k: _drop_image_data(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_drop_image_data(v) for v in obj]
    return obj


def save_log(messages, label):
    """Full messages array as JSON, written on every run (success or not)."""
    LOG_DIR.mkdir(exist_ok=True)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    path = LOG_DIR / f"{label}_run_{ts}.json"
    # Blocks are pydantic objects; dump them as the API sent them
    dumped = json.loads(json.dumps(messages, default=lambda o: o.model_dump(mode="json")))
    path.write_text(json.dumps(_drop_image_data(dumped), indent=2))
    return path


def _yes(question):
    """y/N prompt on the terminal; anything but 'y' is a no."""
    return input(f"\n{question}\n[y/N] ").strip().lower() == "y"


def _approved(code):
    """--confirm gate: show the code and let the user veto it."""
    return _yes(f"--- run_python wants to run ---\n{code}\n-------------------------------\nRun it?")


def _print_usage(totals):
    cost = sum(totals[k] * PRICE[k] for k in PRICE) / 1e6
    print(f"[usage] in={totals['input']:,} out={totals['output']:,} "
          f"cache_read={totals['cache_read']:,} cache_write={totals['cache_write']:,}  ~${cost:.4f}")


def run_agent(task, label="run", max_iterations=20, system_prompt=None, confirm=False):
    client = anthropic.Anthropic()
    kernel = PythonKernel()   # starts on the first run_python call
    system = system_prompt or load_prompt()
    messages = [{"role": "user", "content": task}]
    totals = dict.fromkeys(PRICE, 0)

    try:
        for i in range(1, max_iterations + 1):
            print(f"\n=== step {i}/{max_iterations} ===")
            # Beta endpoint for server-side fallback if a safety classifier declines the request
            with client.beta.messages.stream(
                model=MODEL,
                max_tokens=MAX_TOKENS,
                system=system,
                tools=TOOLS,
                messages=messages,
                thinking={"type": "adaptive", "display": "summarized"},
                cache_control={"type": "ephemeral"},   # caches tools + system + history prefix
                betas=["server-side-fallback-2026-07-01"],
                fallbacks="default",
            ) as stream:
                for text in stream.text_stream:
                    print(text, end="", flush=True)
                response = stream.get_final_message()
            print()

            u = response.usage
            totals["input"] += u.input_tokens
            totals["output"] += u.output_tokens
            totals["cache_read"] += u.cache_read_input_tokens or 0
            totals["cache_write"] += u.cache_creation_input_tokens or 0

            # Thinking blocks must go back unchanged, so append the content as-is
            messages.append({"role": "assistant", "content": response.content})

            if response.stop_reason == "end_turn":
                break
            if response.stop_reason in ("max_tokens", "refusal"):
                print(f"[stopped] stop_reason={response.stop_reason}")
                break

            results = []
            for block in response.content:
                if block.type != "tool_use":
                    continue
                print(f"[tool] {block.name}({json.dumps(block.input)})")
                try:
                    if confirm and block.name == "run_python" and not _approved(block.input["code"]):
                        raise RuntimeError("the user declined to run this code")
                    results.append({"type": "tool_result", "tool_use_id": block.id,
                                    "content": run_tool(block.name, block.input, kernel, ask=_yes)})
                except Exception as e:   # report to the model so it can adjust
                    results.append({"type": "tool_result", "tool_use_id": block.id,
                                    "content": f"Error: {e}", "is_error": True})
                    print(f"[tool error] {e}")
            # All results of one turn go back in a single user message
            messages.append({"role": "user", "content": results})
        else:
            print(f"[stopped] hit max_iterations={max_iterations}")
    finally:
        kernel.shutdown()
        _print_usage(totals)
        print(f"[log] {save_log(messages, label)}")

    return messages
