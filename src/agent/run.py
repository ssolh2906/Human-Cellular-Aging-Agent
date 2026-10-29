#ST
"""
Entry point for the agent.

Usage (from the project root, API key in .env):
  python src/agent/run.py "describe data/Sample_3444.h5ad"
  python src/agent/run.py --sample 3444 --max-iterations 10 "describe this sample"
  python src/agent/run.py --confirm "plot counts per cell for data/Sample_3444.h5ad"
"""

# IMPORTS
import argparse
import os
import sys

from dotenv import load_dotenv

from loop import PROJECT_ROOT, run_agent


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("task", help="what the agent should do")
    parser.add_argument("--sample", help="sample id, only used to name the log file")
    parser.add_argument("--max-iterations", type=int, default=20, help="hard stop on loop steps (default 20)")
    parser.add_argument("--confirm", action="store_true", help="show each run_python code and ask before running it")
    args = parser.parse_args()

    load_dotenv(PROJECT_ROOT / ".env")
    if not os.environ.get("ANTHROPIC_API_KEY"):
        sys.exit("ANTHROPIC_API_KEY not set. Put it in .env (see README).")

    label = f"Sample_{args.sample}" if args.sample else "run"
    run_agent(args.task, label=label, max_iterations=args.max_iterations, confirm=args.confirm)


if __name__ == "__main__":
    main()
