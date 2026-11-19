"""Regenerates every table and figure of the study from a clean checkout."""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEPS = ["primary_results.py", "make_primary_figures.py", "application_example.py", "rheometry_workflow.py"]

for step in STEPS:
    print("running", step, flush=True)
    subprocess.run([sys.executable, str(HERE / step)], cwd=HERE, check=True)
