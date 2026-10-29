#ST
"""
Persistent Python kernel for run_python: variables (e.g. `adata`) survive between calls,
so the agent doesn't reload data every step. Started lazily on first use.
"""

# IMPORTS
import os
import queue
import re
import time
from pathlib import Path

from jupyter_client import KernelManager

PROJECT_ROOT = Path(__file__).resolve().parents[2]
HALF = 6000   # keep this many chars from the head and the tail of long output
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def _clip(text):
    """Head + tail: the end of a traceback or a printout is usually the useful part."""
    if len(text) <= 2 * HALF:
        return text
    return f"{text[:HALF]}\n... [{len(text) - 2 * HALF} chars cut] ...\n{text[-HALF:]}"


class PythonKernel:
    def __init__(self):
        self.km = None
        self.kc = None

    def start(self):
        self.km = KernelManager(kernel_name="python3")   # the kernel installed in this venv
        # Agg: figures are saved to files, never shown in a window
        self.km.start_kernel(cwd=str(PROJECT_ROOT), env={**os.environ, "MPLBACKEND": "Agg"})
        self.kc = self.km.client()
        self.kc.start_channels()
        self.kc.wait_for_ready(timeout=60)

    def shutdown(self):
        if self.km is None:
            return
        self.kc.stop_channels()
        self.km.shutdown_kernel(now=True)
        self.km = self.kc = None

    def execute(self, code, timeout=120):
        """Runs code; returns stdout/stderr, the last expression value and any traceback as text."""
        if self.km is None:
            self.start()
        msg_id = self.kc.execute(code)
        out = []
        deadline = time.time() + timeout
        interrupted = False

        while True:
            remaining = deadline - time.time()
            if remaining <= 0:
                if interrupted:
                    # Interrupt didn't work (e.g. stuck in C code): restart, state is lost
                    self.km.restart_kernel(now=True)
                    self.kc.wait_for_ready(timeout=60)
                    out.append("[kernel did not stop after interrupt and was restarted: all variables are lost]")
                    break
                self.km.interrupt_kernel()
                interrupted = True
                deadline = time.time() + 10   # grace period to reach idle
                out.append(f"[timeout after {timeout}s: execution interrupted, variables are kept]\n")
                continue
            try:
                msg = self.kc.get_iopub_msg(timeout=remaining)
            except queue.Empty:
                continue
            if msg["parent_header"].get("msg_id") != msg_id:
                continue   # leftover output from an earlier call
            kind, content = msg["msg_type"], msg["content"]
            if kind == "stream":
                out.append(content["text"])
            elif kind in ("execute_result", "display_data"):
                out.append(content["data"].get("text/plain", "") + "\n")
            elif kind == "error":
                out.append(ANSI.sub("", "\n".join(content["traceback"])))
            elif kind == "status" and content["execution_state"] == "idle":
                break

        return _clip("".join(out).strip()) or "(no output)"
