"""
Tests for gui/runner.py

Uses real subprocesses (short Python one-liners) — no mocking needed.
Run: pytest tests/gui/test_runner.py -v
"""

import queue
import sys
import threading

import pytest
from gui.runner import run_pipeline, start_pipeline_thread


def _collect(cmd) -> list:
    """Run run_pipeline synchronously and return all queued messages."""
    q = queue.Queue()
    run_pipeline(cmd, q, [None])
    msgs = []
    while not q.empty():
        msgs.append(q.get())
    return msgs


class TestRunPipeline:
    def test_successful_command_exit_zero(self):
        msgs = _collect([sys.executable, "-c", "pass"])
        exits = [m for m in msgs if m[0] == "exit"]
        assert exits == [("exit", 0)]

    def test_output_lines_forwarded(self):
        msgs = _collect([sys.executable, "-c", "print('hello')"])
        lines = [m[1] for m in msgs if m[0] == "line"]
        assert "hello" in lines

    def test_multiple_lines_in_order(self):
        msgs = _collect([sys.executable, "-c",
                         "print('one'); print('two'); print('three')"])
        lines = [m[1] for m in msgs if m[0] == "line"]
        assert lines == ["one", "two", "three"]

    def test_failing_command_nonzero_exit(self):
        msgs = _collect([sys.executable, "-c", "raise SystemExit(3)"])
        exits = [m for m in msgs if m[0] == "exit"]
        assert exits[0] == ("exit", 3)

    def test_exit_is_last_message(self):
        msgs = _collect([sys.executable, "-c",
                         "print('a'); print('b')"])
        assert msgs[-1][0] == "exit"

    def test_stderr_merged_into_lines(self):
        """stderr is merged into stdout so errors appear in the log."""
        msgs = _collect([sys.executable, "-c",
                         "import sys; sys.stderr.write('err_line\\n')"])
        lines = [m[1] for m in msgs if m[0] == "line"]
        assert any("err_line" in l for l in lines)

    def test_invalid_binary_sends_error_and_exit(self):
        msgs = _collect(["/nonexistent_binary_vi_pipeline"])
        types = [m[0] for m in msgs]
        assert "error" in types
        assert "exit" in types
        exit_code = next(m[1] for m in msgs if m[0] == "exit")
        assert exit_code != 0

    def test_no_output_command_has_only_exit(self):
        msgs = _collect([sys.executable, "-c", "pass"])
        non_exit = [m for m in msgs if m[0] != "exit"]
        assert non_exit == []

    def test_exit_code_preserved_for_code_1(self):
        msgs = _collect([sys.executable, "-c", "raise SystemExit(1)"])
        exits = [m for m in msgs if m[0] == "exit"]
        assert exits[0] == ("exit", 1)


class TestStartPipelineThread:
    def test_returns_thread_instance(self):
        q = queue.Queue()
        t = start_pipeline_thread([sys.executable, "-c", "pass"], q, [None])
        assert isinstance(t, threading.Thread)
        t.join(timeout=5)

    def test_thread_is_daemon(self):
        q = queue.Queue()
        t = start_pipeline_thread([sys.executable, "-c", "pass"], q, [None])
        assert t.daemon is True
        t.join(timeout=5)

    def test_thread_completes(self):
        q = queue.Queue()
        t = start_pipeline_thread([sys.executable, "-c", "pass"], q, [None])
        t.join(timeout=5)
        assert not t.is_alive()

    def test_thread_posts_exit_to_queue(self):
        q = queue.Queue()
        t = start_pipeline_thread([sys.executable, "-c", "pass"], q, [None])
        t.join(timeout=5)
        msgs = []
        while not q.empty():
            msgs.append(q.get())
        exits = [m for m in msgs if m[0] == "exit"]
        assert len(exits) == 1

    def test_proc_holder_populated(self):
        """proc_holder[0] is set to the Popen object once the process starts."""
        import subprocess
        q = queue.Queue()
        holder = [None]
        t = start_pipeline_thread(
            [sys.executable, "-c", "import time; time.sleep(0.2)"], q, holder
        )
        t.join(timeout=5)
        assert isinstance(holder[0], subprocess.Popen)

    def test_stop_via_proc_holder(self):
        """Terminating proc_holder[0] causes the thread to exit and post exit code."""
        import time
        q = queue.Queue()
        holder = [None]
        t = start_pipeline_thread(
            [sys.executable, "-c", "import time; time.sleep(60)"], q, holder
        )
        # Wait until the process actually starts
        deadline = time.monotonic() + 5
        while holder[0] is None and time.monotonic() < deadline:
            time.sleep(0.05)
        assert holder[0] is not None, "Process never started"
        holder[0].terminate()
        t.join(timeout=5)
        msgs = []
        while not q.empty():
            msgs.append(q.get())
        exits = [m for m in msgs if m[0] == "exit"]
        assert len(exits) == 1
        assert exits[0][1] != 0
