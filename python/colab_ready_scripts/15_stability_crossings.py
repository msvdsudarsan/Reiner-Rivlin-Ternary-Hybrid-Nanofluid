# Stability crossings at K=+0.3 and K=-0.3 (S=1.5), evidence table and figure.
# Needs the folder ../rr_constitutive (stability_tracked.py etc.); run from this folder.
import os, sys, runpy
here = os.path.dirname(os.path.abspath(__file__)); rc = os.path.join(here, "..", "rr_constitutive")
sys.path.insert(0, rc); os.chdir(rc)
for s in ("run_stability_check.py", "mode_evidence.py", "make_stability_figure.py"):
    runpy.run_path(os.path.join(rc, s), run_name="__main__")
