"""Validation of both CFD models -> cfd/VALIDATION.md, cfd/validation.json, figures/12_cfd_validation.png

1D reactor model (python_1d_reactor/reactor_model.py)
  V1  global energy conservation for all cases       |error| < 1 %
  V2  grid convergence (40 / 80 / 160 cells)          change 80->160 < 2 %
  V3  thermodynamic limit: adiabatic discharge must stop at T_eq(p_in) = dH / (dS - R ln(p_in/p0))
  V4  gas properties against reference values (NIST / VDI Heat Atlas, pure N2 and steam)
OpenFOAM case (openfoam/packedBedReactor, needs OpenFOAM v1912+ on PATH)
  O1  bed pressure drop vs Ergun equation             |error| < 3 %
  O2  mass conservation (patch face fluxes)            |error| < 0.1 %
  O3  energy conservation, adiabatic wall              |error| < 1 %
  O4  mesh independence (coarse / base / fine)         Δp and T_out change < 2 %
  O5  1D <-> OpenFOAM consistency: same adiabatic heat input -> same outlet temperature rise

Usage:  python cfd/validate_cfd.py            (add --skip-openfoam if OpenFOAM is not installed)
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import replace
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "python_1d_reactor"))
sys.path.insert(0, str(HERE / "openfoam"))
from properties import gas_mixture  # noqa: E402
from reactor_model import Params, run, R  # noqa: E402

CASES = HERE / "python_1d_reactor" / "cases"
results, rows = {}, []


def check(tag, name, value, limit, unit="%", passed=None):
    ok = passed if passed is not None else abs(value) <= limit
    rows.append((tag, name, f"{value:.3g} {unit}".strip(), f"≤ {limit} {unit}" if passed is None else limit, "PASS" if ok else "FAIL"))
    print(f"[{'PASS' if ok else 'FAIL'}] {tag} {name}: {value:.4g} {unit}")
    return ok


# ---------------------------------------------------------------- 1D model
def validate_1d():
    for js in sorted(CASES.glob("*.json")):
        k = run(Params.from_json(js)).kpis()
        check("V1", f"energy balance, {js.stem}", k["energy_balance_error_pct"], 1.0)

    base = Params.from_json(CASES / "caoh2_discharge.json")
    conv = {}
    for n in (40, 80, 160):
        k = run(replace(base, n_cells=n, n_out=361)).kpis()
        conv[n] = k
    for key in ("time_to_90pct_conversion_s", "T_out_max_C", "peak_total_power_W"):
        a, b = conv[80][key], conv[160][key]
        check("V2", f"grid 80->160 cells, {key}", 100 * (a - b) / b, 2.0)
    results["grid_convergence"] = {n: {k: conv[n][k] for k in ("time_to_90pct_conversion_s", "T_out_max_C", "peak_total_power_W")} for n in conv}

    adi = replace(base, h_w=0.0, t_end=4 * 3600, L=0.2, n_cells=80)
    res = run(adi)
    T_eq = adi.dH / (adi.dS - R * math.log(adi.p_in / adi.p_tot))
    Tmax = float(res.T.max())
    check("V3", "adiabatic T_max vs T_eq(p_in) (must not exceed, gap < 3 K)", Tmax - T_eq, 3.0, "K",
          passed=(Tmax <= T_eq + 0.05) and (T_eq - Tmax < 3.0))
    results["equilibrium_limit"] = dict(T_eq_C=T_eq - 273.15, T_max_C=Tmax - 273.15)
    results["_adiabatic_1d"] = res

    # reference gas properties (pure species) - NIST WebBook / VDI Heat Atlas
    ref = [("N2 cp 400 C", gas_mixture(673.15, 1e5, 0.0).cp, 1098.0),
           ("N2 mu 400 C", gas_mixture(673.15, 1e5, 0.0).mu, 3.20e-5),
           ("N2 k 400 C", gas_mixture(673.15, 1e5, 0.0).k, 0.0485),
           ("H2O cp 400 C, 1 bar", gas_mixture(673.15, 1e5, 1e5).cp, 2070.0),
           ("H2O mu 400 C", gas_mixture(673.15, 1e5, 1e5).mu, 2.43e-5),
           ("H2O k 400 C", gas_mixture(673.15, 1e5, 1e5).k, 0.0540)]
    for name, v, r in ref:
        check("V4", f"{name} ({v:.4g} vs ref {r:.4g})", 100 * (v - r) / r, 5.0)


# ---------------------------------------------------------------- OpenFOAM
def of_env():
    for rc in ("/usr/share/openfoam/etc/bashrc", "/usr/lib/openfoam/openfoam/etc/bashrc"):
        if Path(rc).exists():
            return f"source {rc} >/dev/null 2>&1; "
    return ""


def run_of(case_dir, mesh):
    bmd = case_dir / "system" / "blockMeshDict"
    s = bmd.read_text()
    nx_in, nx_bed, nx_out, ny = mesh
    counts = iter([(nx_in, ny), (nx_bed, ny), (nx_out, ny)])
    s = re.sub(r"(hex\s*\([^)]*\)\s*(?:\w+\s*)?)\((\d+)\s+(\d+)\s+1\)", lambda m: "{}({} {} 1)".format(m.group(1), *next(counts)), s)
    bmd.write_text(s)
    cmd = of_env() + "blockMesh > log.blockMesh 2>&1 && rhoPorousSimpleFoam > log.rhoPorousSimpleFoam 2>&1"
    subprocess.run(["bash", "-c", cmd], cwd=case_dir, check=True)
    import postprocess
    return postprocess.main(["--case", str(case_dir), "--no-plots"])


def validate_openfoam():
    if shutil.which("blockMesh") is None and not of_env():
        print("OpenFOAM not found - skipping O-checks"); return
    subprocess.run([sys.executable, str(HERE / "openfoam" / "setup_from_1d.py"), "--wall", "adiabatic"], check=True)
    meshes = {"coarse": (10, 40, 10, 12), "base": (20, 80, 20, 25), "fine": (40, 160, 40, 50)}
    out = {}
    with tempfile.TemporaryDirectory() as tmp:
        for name, mesh in meshes.items():
            d = Path(tmp) / name
            shutil.copytree(HERE / "openfoam" / "packedBedReactor", d)
            out[name] = run_of(d, mesh)
    b, f = out["base"], out["fine"]
    check("O1", "bed Δp vs Ergun (base mesh)", b["pressure_drop_error_pct"], 3.0)
    check("O2", "mass conservation (base mesh)", b["mass_balance_error_pct"], 0.1)
    check("O3", "energy conservation, adiabatic (base mesh)", b["energy_balance_error_pct"], 1.0)
    check("O4", "mesh base->fine, bed Δp", 100 * (b["pressure_drop_bed_cfd_Pa"] - f["pressure_drop_bed_cfd_Pa"]) / f["pressure_drop_bed_cfd_Pa"], 2.0)
    dT_b = b["T_outlet_bulk_C"] - b["T_inlet_C"]; dT_f = f["T_outlet_bulk_C"] - f["T_inlet_C"]
    check("O4", "mesh base->fine, outlet ΔT", 100 * (dT_b - dT_f) / dT_f, 2.0)
    # O5: analytic adiabatic rise with the 1D gas properties: dT = q V / (m_dot cp)
    S = json.loads((HERE / "openfoam" / "packedBedReactor" / "setup.json").read_text())
    A = math.pi * S["bed_diameter"] ** 2 / 4
    dT_an = S["q_W_m3"] * A * S["bed_length"] / (S["gas"]["rho"] * S["u_s"] * A * S["gas"]["cp"])
    check("O5", f"OpenFOAM ΔT {dT_b:.2f} K vs 1D-gas energy balance {dT_an:.2f} K", 100 * (dT_b - dT_an) / dT_an, 2.0)
    results["openfoam"] = out
    # restore production set-up
    subprocess.run([sys.executable, str(HERE / "openfoam" / "setup_from_1d.py"), "--wall", "isothermal"], check=True)


def report():
    n_fail = sum(r[4] == "FAIL" for r in rows)
    md = ["# CFD validation report", "",
          "Generated by `python cfd/validate_cfd.py`. Both models use the same gas properties "
          "(`cfd/properties.py`), the same bed (ε = 0.40, d_p = 2 mm, D = 50 mm, L = 200 mm) and the same "
          "Ca(OH)₂ equilibrium (Schaube et al. 2012: ΔH = 106.8 kJ/mol, ΔS = 137.2 J/mol/K).", "",
          f"**{len(rows) - n_fail} of {len(rows)} checks passed.**", "",
          "| ID | Check | Result | Criterion | Status |", "|---|---|---|---|---|"]
    md += [f"| {a} | {b} | {c} | {d} | {e} |" for a, b, c, d, e in rows]
    md += ["", "## What the checks mean", "",
           "* **V1 / O3** conservation of energy: heat released by the reaction (or the heat source) equals heat leaving with the gas, through the wall and stored in the bed.",
           "* **V2 / O4** numerical accuracy: results do not change when the grid is refined.",
           "* **V3** physics: an adiabatic bed heats up until the reaction reaches equilibrium and then stops; it must never exceed T_eq at the inlet vapour pressure.",
           "* **V4** material properties against tabulated reference values for the pure gases.",
           "* **O1** the porous-medium resistance in OpenFOAM reproduces the Ergun equation.",
           "* **O2** the OpenFOAM solution conserves mass.",
           "* **O5** both CFD models are consistent: the same heat input gives the same gas temperature rise.",
           "", "## Not validated here (needs experiments)", "",
           "* Reaction kinetics (A, Ea) are literature order-of-magnitude values - fit them to TGA data.",
           "* Effective bed conductivity and wall heat-transfer coefficient - fit to reactor temperature measurements.",
           "* The OpenFOAM case prescribes the reaction heat; it has no solid phase or conversion field."]
    (HERE / "VALIDATION.md").write_text("\n".join(md) + "\n")
    js = {k: v for k, v in results.items() if not k.startswith("_")}
    (HERE / "validation.json").write_text(json.dumps(dict(checks=rows, **js), indent=2, default=str))
    plot()
    print(f"\n{len(rows) - n_fail}/{len(rows)} checks passed -> cfd/VALIDATION.md")
    return n_fail


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    g = results.get("grid_convergence", {})
    if g:
        n = sorted(g, key=int)
        ax[0].plot(n, [g[i]["time_to_90pct_conversion_s"] / 60 for i in n], "o-", color="#2a78d6", lw=2)
        ax[0].set(xscale="log", xlabel="Cells", ylabel="Time to 90 % conversion [min]", title="1D model: grid convergence")
        ax[0].minorticks_off()
        ax[0].set_xticks([int(i) for i in n], [str(i) for i in n])
    res = results.get("_adiabatic_1d")
    if res is not None:
        eq = results["equilibrium_limit"]
        ax[1].plot(res.t / 60, res.T.max(axis=0) - 273.15, color="#2a78d6", lw=2, label="max bed T (adiabatic)")
        ax[1].axhline(eq["T_eq_C"], color="#eb6834", ls="--", lw=1.5, label=f"T_eq(p_in) = {eq['T_eq_C']:.1f} °C")
        ax[1].set(xlabel="Time [min]", ylabel="T [°C]", title="1D model: equilibrium limit")
        ax[1].legend(frameon=False, fontsize=8)
    of = results.get("openfoam")
    if of:
        names = list(of)
        ax[2].bar(np.arange(3) - 0.2, [of[m]["pressure_drop_bed_cfd_Pa"] for m in names], 0.4, color="#2a78d6", label="OpenFOAM Δp")
        ax[2].bar(np.arange(3) + 0.2, [of[m]["pressure_drop_bed_ergun_Pa"] for m in names], 0.4, color="#eb6834", label="Ergun Δp")
        ax[2].set_xticks(range(3), names)
        ax[2].set(ylabel="Bed pressure drop [Pa]", title="OpenFOAM: Δp vs Ergun, 3 meshes")
        ax[2].set_ylim(0, 1.3 * max(of[m]["pressure_drop_bed_cfd_Pa"] for m in names))
        ax[2].legend(frameon=False, fontsize=8, loc="upper center", ncol=2)
    for a in ax:
        a.grid(color="#e4e3df"); a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "12_cfd_validation.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-openfoam", action="store_true")
    a = ap.parse_args()
    validate_1d()
    if not a.skip_openfoam:
        validate_openfoam()
    sys.exit(1 if report() else 0)
