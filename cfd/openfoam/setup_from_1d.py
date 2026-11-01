"""Write consistent material properties and operating conditions into the OpenFOAM case.

Everything comes from the same sources as the 1D model, so both CFD models describe the
same reactor:
  * gas mixture (humid N2) M, cp, mu, Pr  -> constant/thermophysicalProperties   (cfd/properties.py)
  * Ergun -> Darcy-Forchheimer d, f from the 1D case's eps and dp -> constant/porosityProperties
  * volumetric reaction heat from the 1D results                 -> constant/fvOptions
  * inlet velocity, inlet / wall temperature, outlet pressure    -> 0/U, 0/T, 0/p
  * wall: 'isothermal' (heat-transfer fluid at T_wall) or 'adiabatic' (validation: closed energy balance)

Usage:
    python setup_from_1d.py                                   # Ca(OH)2 discharge case, isothermal wall
    python setup_from_1d.py --wall adiabatic                  # validation set-up
    python setup_from_1d.py --case-json ../python_1d_reactor/cases/caoh2_charge.json
"""
import argparse
import json
import math
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent / "python_1d_reactor"))
from properties import gas_mixture  # noqa: E402
from reactor_model import Params, run  # noqa: E402

CASE = HERE / "packedBedReactor"


def sub(path, pattern, repl, count=0):
    s = path.read_text()
    s2, n = re.subn(pattern, repl, s, count=count, flags=re.S)
    if n == 0:
        raise SystemExit(f"pattern not found in {path}: {pattern}")
    path.write_text(s2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case-json", default=str(HERE.parent / "python_1d_reactor" / "cases" / "caoh2_discharge.json"))
    ap.add_argument("--wall", choices=["isothermal", "adiabatic"], default="isothermal")
    ap.add_argument("--source", choices=["gas", "total"], default="gas",
                    help="gas: heat the 1D model removes with the gas (default, because this OpenFOAM model has no "
                         "solid-phase conduction to the wall); total: all reaction heat")
    ap.add_argument("--q", type=float, help="override volumetric heat release [W/m3]")
    a = ap.parse_args()

    P = Params.from_json(a.case_json)
    res = run(P)                      # run the 1D model so the heat source is always up to date
    k = res.kpis()
    g = gas_mixture(P.T_in, P.p_tot, P.p_in)

    # gas properties
    tp = CASE / "constant" / "thermophysicalProperties"
    sub(tp, r"molWeight\s+[0-9.eE+-]+;", f"molWeight {g.M * 1e3:.4f};")
    sub(tp, r"Cp\s+[0-9.eE+-]+;", f"Cp {g.cp:.2f};")
    sub(tp, r"mu\s+[0-9.eE+-]+;", f"mu {g.mu:.5e};")
    sub(tp, r"Pr\s+[0-9.eE+-]+;", f"Pr {g.Pr:.4f};")
    sub(tp, r"// GAS:.*?\n", f"// GAS: y_H2O = {g.y_H2O:.2f}, T = {P.T_in - 273.15:.0f} C (cfd/properties.py)\n")

    # Ergun -> Darcy-Forchheimer
    d = 150 * (1 - P.eps) ** 2 / (P.eps ** 3 * P.dp ** 2)
    f = 3.5 * (1 - P.eps) / (P.eps ** 3 * P.dp)
    pp = CASE / "constant" / "porosityProperties"
    vec = lambda v: f"({v:.6e} {v:.6e} {v:.6e})"
    sub(pp, r"\bd\s+\([^)]*\);", f"d   {vec(d)};")
    sub(pp, r"\bf\s+\([^)]*\);", f"f   {vec(f)};")
    sub(pp, r"dp = [^,]*, eps = [0-9.]+", f"dp = {P.dp * 1000:g} mm, eps = {P.eps:.2f}")

    # heat source
    vol = math.pi * P.D ** 2 / 4 * P.L
    q = a.q if a.q is not None else k["mean_gas_power_W" if a.source == "gas" else "mean_total_power_W"] / vol
    if P.mode == "charge":
        q = -abs(q)
    sub(CASE / "constant" / "fvOptions", r"\bh\s+\([^)]*\);", f"h   ({q:.4e} 0);")

    # boundary / initial conditions
    sub(CASE / "0" / "U", r"uniform \([0-9.eE+-]+ 0 0\)", f"uniform ({P.u_s} 0 0)")
    sub(CASE / "0" / "T", r"internalField\s+uniform [0-9.]+;", f"internalField   uniform {P.T_in};")
    sub(CASE / "0" / "T", r"(inlet\s*\{[^}]*?)uniform [0-9.]+", rf"\g<1>uniform {P.T_in}")
    sub(CASE / "0" / "T", r"inletValue uniform [0-9.]+; value uniform [0-9.]+",
        f"inletValue uniform {P.T_in}; value uniform {P.T_in}")
    wall = (f"wall         {{ type fixedValue; value uniform {P.T_wall}; }}   // heat-transfer-fluid side"
            if a.wall == "isothermal" else "wall         { type zeroGradient; }   // adiabatic (validation)")
    sub(CASE / "0" / "T", r"wall\s*\{[^}]*\}[^\n]*", wall)
    for fld in ("p",):
        sub(CASE / "0" / fld, r"internalField\s+uniform [0-9.eE+-]+;", f"internalField   uniform {P.p_tot:g};")
        sub(CASE / "0" / fld, r"(outlet\s*\{[^}]*?)uniform [0-9.eE+-]+", rf"\g<1>uniform {P.p_tot:g}")

    setup = dict(case_json=a.case_json, wall=a.wall, source=a.source, q_W_m3=q, eps=P.eps, dp=P.dp,
                 darcy_d=d, forchheimer_f=f, u_s=P.u_s, T_in=P.T_in, T_wall=P.T_wall, p_tot=P.p_tot,
                 bed_length=P.L, bed_diameter=P.D, gas=dict(M=g.M, rho=g.rho, cp=g.cp, mu=g.mu, k=g.k, Pr=g.Pr))
    (CASE / "setup.json").write_text(json.dumps(setup, indent=2))
    print(f"gas: M={g.M * 1e3:.2f} g/mol cp={g.cp:.0f} J/kgK mu={g.mu:.3e} Pa s Pr={g.Pr:.3f} | "
          f"Ergun d={d:.3e} f={f:.3e} | q={q:.3e} W/m3 | wall={a.wall}")


if __name__ == "__main__":
    main()
