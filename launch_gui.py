#!/usr/bin/env python3
"""
launch_gui.py — Start the VI_Pipeline desktop GUI.

Usage:
    python launch_gui.py
"""

import sys
from pathlib import Path

# Ensure the project root is on sys.path when launched from any directory.
sys.path.insert(0, str(Path(__file__).parent))

from gui.app import App

if __name__ == "__main__":
    App().mainloop()
