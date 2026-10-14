"""Rebuild from archived data, audit, then render. Run with Python 3.12."""

from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
for script in [
    ROOT / "Fig4_panel_c" / "generate_fig4_panel_c.py",
    ROOT / "code" / "build_figure.py",
    ROOT / "code" / "verify_evidence.py",
    ROOT / "code" / "render_delivery.py",
]:
    subprocess.run([sys.executable, str(script)], check=True)
