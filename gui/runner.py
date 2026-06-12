"""
runner.py — Run the VI_Pipeline CLI in a background thread.

Messages posted to the queue:
  ("line",  str)  — one line of stdout/stderr output
  ("exit",  int)  — process return code (always the last message)
  ("error", str)  — internal error before the process could start
"""

import os
import queue
import subprocess
import threading


def run_pipeline(cmd: list, q: queue.Queue, proc_holder: list) -> None:
    """
    Run *cmd* in a subprocess, forwarding all output to *q*.

    proc_holder is a one-element list [None] that is populated with the
    subprocess.Popen object as soon as it starts.
    """
    try:
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,   # merge stderr so the log is in order
            text=True,
            bufsize=1,
            preexec_fn=os.setsid,       # new process group — lets us kill the whole tree
        )
        proc_holder[0] = proc
        for line in proc.stdout:
            q.put(("line", line.rstrip("\n")))
        proc.wait()
        q.put(("exit", proc.returncode))
    except Exception as exc:           # noqa: BLE001
        q.put(("error", str(exc)))
        q.put(("exit", 1))


def start_pipeline_thread(
    cmd: list,
    q: queue.Queue,
    proc_holder: list,
) -> threading.Thread:
    """
    Start *run_pipeline* in a daemon thread and return it.

    proc_holder must be a one-element list initialised to [None].
    It will be populated with the Popen object once the process starts,
    allowing the caller to terminate it via proc_holder[0].terminate().
    """
    t = threading.Thread(target=run_pipeline, args=(cmd, q, proc_holder), daemon=True)
    t.start()
    return t
