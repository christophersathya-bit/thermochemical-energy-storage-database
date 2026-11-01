"""Rebuild everything: properties -> SQLite database -> classification -> figures -> ML -> 1D reactor cases."""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    ["scripts/compute_properties.py"],
    ["database/build_database.py"],
    ["scripts/verification_report.py"],
    ["scripts/build_thermophysical.py"],
    ["scripts/classification.py"],
    ["scripts/plots.py"],
    ["cfd/python_1d_reactor/reactor_model.py", "--config", "cfd/python_1d_reactor/cases/caoh2_discharge.json"],
    ["cfd/python_1d_reactor/reactor_model.py", "--config", "cfd/python_1d_reactor/cases/caoh2_charge.json"],
    ["cfd/python_1d_reactor/reactor_model.py", "--config", "cfd/python_1d_reactor/cases/mgso4_discharge.json"],
    ["ml/tces_ml.py", "all"],
]
if "--validate" in sys.argv:
    STEPS.append(["cfd/validate_cfd.py"] + ([] if "--with-openfoam" in sys.argv else ["--skip-openfoam"]))
if "--with-surrogate" in sys.argv:
    STEPS.append(["ml/reactor_surrogate.py", "--n", "300"])

for s in STEPS:
    print(f"\n=== {' '.join(s)}")
    subprocess.run([sys.executable, str(ROOT / s[0]), *s[1:]], check=True, cwd=ROOT)
print("\nAll steps finished.")
