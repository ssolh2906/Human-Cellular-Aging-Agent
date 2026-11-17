#ST
"""
Entry point for the agent.

Usage (from the project root, API key in .env):
  python src/agent/run.py "describe data/Sample_3444.h5ad"
  python src/agent/run.py --sample 3444 --max-iterations 10 "describe this sample"
  python src/agent/run.py --confirm "plot counts per cell for data/Sample_3444.h5ad"

Stages: the task comes from src/agent/prompts/<stage>.md; text given after it is added as extra instructions.
  python src/agent/run.py --stage bubbletree --max-iterations 60
  python src/agent/run.py --stage bubbletree "Only process samples 3444 and 1931 for now."
"""

# IMPORTS
import argparse
import os
import sys

from dotenv import load_dotenv

from loop import PROJECT_ROOT, PROMPT_DIR, load_prompt, run_agent


def main():
    # Windows defaults to cp1252 when output is piped; one non-Latin character from the model would end the run
    sys.stdout.reconfigure(encoding="utf-8")
    # Every prompts/*.md except the shared system prompt is a stage, so a new file needs no code change
    stages = sorted(p.stem for p in PROMPT_DIR.glob("*.md") if p.stem != "system")

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("task", nargs="?", help="what the agent should do (with --stage: extra instructions)")
    parser.add_argument("--stage", choices=stages, help="take the task from prompts/<stage>.md")
    parser.add_argument("--sample", help="sample id, only used to name the log file")
    parser.add_argument("--max-iterations", type=int, default=20, help="hard stop on loop steps (default 20)")
    parser.add_argument("--confirm", action="store_true", help="show each run_python code and ask before running it")
    args = parser.parse_args()
    if not (args.task or args.stage):
        parser.error("give a task, a --stage, or both")

    load_dotenv(PROJECT_ROOT / ".env")
    if not os.environ.get("ANTHROPIC_API_KEY"):
        sys.exit("ANTHROPIC_API_KEY not set. Put it in .env (see README).")

    stage_prompt = load_prompt(args.stage) if args.stage else None
    task = "\n\n".join(part for part in (stage_prompt, args.task) if part)
    label = f"Sample_{args.sample}" if args.sample else args.stage or "run"
    run_agent(task, label=label, max_iterations=args.max_iterations, confirm=args.confirm)


if __name__ == "__main__":
    main()
