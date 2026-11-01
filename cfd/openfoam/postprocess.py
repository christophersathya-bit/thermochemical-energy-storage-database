"""Post-process and validate the OpenFOAM packed-bed case (axisymmetric wedge).

Reads p, T, U from the latest time directory, rebuilds the blockMesh cell layout
(inlet plenum | bed | outlet plenum, x-fastest then radial) and checks:
  1. pressure drop over the bed        vs Ergun equation (same eps, dp, mu, rho)
  2. mass conservation                 inlet vs outlet mass flow
  3. energy conservation (adiabatic)   q * V_bed  vs  m_dot * cp * (T_out,bulk - T_in)
     (isothermal wall: the difference is reported as heat removed through the wall)
Writes packedBedReactor/postProcessing_py/{summary.json, fields_2d.png, axial_profiles.png}

Usage:  python postprocess.py [--case packedBedReactor]
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
X_EDGES = (0.00, 0.05, 0.25, 0.30)                # inlet plenum | bed | outlet plenum [m]
BLOCKS = [(0.00, 0.05, 20), (0.05, 0.25, 80), (0.25, 0.30, 20)]  # overwritten from blockMeshDict
NY = 25


def read_blocks(case):
    """Cell counts per block from system/blockMeshDict."""
    global BLOCKS, NY
    s = (case / "system" / "blockMeshDict").read_text()
    n = re.findall(r"hex\s*\([^)]*\)\s*(?:\w+\s*)?\((\d+)\s+(\d+)\s+1\)", s)
    BLOCKS = [(X_EDGES[i], X_EDGES[i + 1], int(nx)) for i, (nx, _) in enumerate(n)]
    NY = int(n[0][1])
R_UNIV = 8.314462618


def read_internal(path, vector=False):
    s = path.read_text()
    m = re.search(r"internalField\s+nonuniform\s+List<\w+>\s*\n?(\d+)\s*\n?\(", s)
    if not m:
        raise SystemExit(f"{path}: field is uniform - has the solver run?")
    n = int(m.group(1))
    body = s[m.end():]
    if vector:
        vals = re.findall(r"\(([^()]+)\)", body)[:n]
        return np.array([[float(v) for v in t.split()] for t in vals])
    return np.array(body.split(")")[0].split(), dtype=float)[:n]


def read_patch(path, patch, n=None):
    """Values of a patch in boundaryField (uniform or nonuniform scalar list)."""
    s = path.read_text()
    s = s[s.index("boundaryField"):]
    m = re.search(rf"\b{patch}\s*\{{(.*?)\}}", s, re.S)
    blk = m.group(1)
    mu = re.search(r"value\s+uniform\s+([0-9.eE+-]+)", blk)
    if mu:
        return np.full(n or NY, float(mu.group(1)))
    ml = re.search(r"value\s+nonuniform\s+List<scalar>\s*(\d+)\s*\(([^)]*)\)", blk, re.S)
    return np.array(ml.group(2).split(), dtype=float)


def to_grid(field):
    parts, i = [], 0
    for _, _, nx in BLOCKS:
        parts.append(field[i:i + nx * NY].reshape(NY, nx))
        i += nx * NY
    return np.hstack(parts)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=str(HERE / "packedBedReactor"))
    ap.add_argument("--no-plots", action="store_true")
    a = ap.parse_args(argv)
    case = Path(a.case)
    read_blocks(case)
    S = json.loads((case / "setup.json").read_text())
    gas = S["gas"]
    Rr = S["bed_diameter"] / 2
    times = sorted((d for d in case.iterdir() if d.is_dir() and re.fullmatch(r"[0-9.e+-]+", d.name) and float(d.name) > 0),
                   key=lambda d: float(d.name))
    if not times:
        raise SystemExit("No result time directories - run ./Allrun first.")
    t = times[-1]
    p, T = to_grid(read_internal(t / "p")), to_grid(read_internal(t / "T"))
    Ux = to_grid(read_internal(t / "U", vector=True)[:, 0])
    x = np.concatenate([x0 + (np.arange(nx) + 0.5) * (x1 - x0) / nx for x0, x1, nx in BLOCKS])
    dr = Rr / NY
    r = (np.arange(NY) + 0.5) * dr
    dA = 2 * math.pi * r * dr                      # annulus area of each radial cell (full 360 deg)
    A = math.pi * Rr ** 2
    rho = p * (gas["M"]) / (R_UNIV * T)
    mflux = rho * Ux * dA[:, None]                 # kg/s through each annulus

    def area_avg(F, i):
        return float((F[:, i] * dA).sum() / dA.sum())

    def bulk_T(i):
        return float((mflux[:, i] * T[:, i]).sum() / mflux[:, i].sum())

    i_in, i_out = np.searchsorted(x, 0.05), np.searchsorted(x, 0.25) - 1
    # 1. Ergun (evaluated with the mean gas density and superficial velocity in the bed)
    L_c = x[i_out] - x[i_in]                       # distance between first and last bed cell centres
    dp_cfd = area_avg(p, i_in) - area_avg(p, i_out)
    xs = slice(i_in, i_out + 1)
    rho_b = float((rho[:, xs] * dA[:, None]).sum() / dA.sum() / (i_out - i_in + 1))
    u_b = float(mflux[:, xs].sum(axis=0).mean()) / (rho_b * A)
    dp_ergun = (gas["mu"] * S["darcy_d"] * u_b + 0.5 * rho_b * S["forchheimer_f"] * u_b ** 2) * L_c
    # 2. mass balance from the face fluxes phi [kg/s per wedge]
    wedge = 360.0 / 5.0
    phi_in, phi_out = read_patch(t / "phi", "inlet"), read_patch(t / "phi", "outlet")
    m_dot, m_out = float(-phi_in.sum() * wedge), float(phi_out.sum() * wedge)
    # 3. energy balance
    V_bed = A * S["bed_length"]
    Q_src = S["q_W_m3"] * V_bed
    T_in_f, T_out_f = read_patch(t / "T", "inlet"), read_patch(t / "T", "outlet")
    T_in_b = float((phi_in * T_in_f).sum() / phi_in.sum())
    T_out_b = float((phi_out * T_out_f).sum() / phi_out.sum())
    Q_gas = gas["cp"] * (m_out * T_out_b - m_dot * T_in_b)
    summary = dict(
        time=t.name, wall=S["wall"],
        pressure_drop_bed_cfd_Pa=dp_cfd, pressure_drop_bed_ergun_Pa=dp_ergun,
        pressure_drop_error_pct=100 * (dp_cfd - dp_ergun) / dp_ergun,
        mass_flow_in_kg_s=m_dot, mass_flow_out_kg_s=m_out, mass_balance_error_pct=100 * (m_out - m_dot) / m_dot,
        heat_source_W=Q_src, heat_to_gas_W=Q_gas, heat_through_wall_W=Q_src - Q_gas,
        energy_balance_error_pct=100 * (Q_gas - Q_src) / Q_src if S["wall"] == "adiabatic" else None,
        T_inlet_C=T_in_b - 273.15, T_outlet_bulk_C=T_out_b - 273.15, T_max_C=float(T.max() - 273.15),
    )
    out = case / "postProcessing_py"
    out.mkdir(exist_ok=True)
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    if a.no_plots:
        return summary

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(10, 5.5), sharex=True)
    for a_, F, lab, cm in [(ax[0], T - 273.15, "T [°C]", "inferno"), (ax[1], p - S["p_tot"], "p - p_out [Pa]", "viridis")]:
        im = a_.pcolormesh(x * 100, r * 100, F, cmap=cm, shading="auto")
        fig.colorbar(im, ax=a_, label=lab)
        for xb in (5, 25):
            a_.axvline(xb, color="w", ls="--", lw=1)
        a_.set_ylabel("r [cm]")
    ax[0].set_title(f"OpenFOAM rhoPorousSimpleFoam, axisymmetric bed (dashed), wall: {S['wall']}", loc="left", fontweight="bold")
    ax[1].set_xlabel("x [cm]")
    fig.tight_layout(); fig.savefig(out / "fields_2d.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    Tb = (mflux * T).sum(axis=0) / mflux.sum(axis=0)
    ax[0].plot(x * 100, Tb - 273.15, color="#2a78d6", lw=2)
    ax[0].set(xlabel="x [cm]", ylabel="bulk T [°C]", title="Axial bulk temperature")
    ax[1].plot(x * 100, (p * dA[:, None]).sum(axis=0) / dA.sum() - S["p_tot"], color="#2a78d6", lw=2, label="OpenFOAM")
    ax[1].set(xlabel="x [cm]", ylabel="p - p_out [Pa]", title=f"Pressure: bed Δp {dp_cfd:.0f} Pa (Ergun {dp_ergun:.0f} Pa)")
    for a_ in ax:
        a_.axvspan(5, 25, color="#e4e3df", alpha=.5, lw=0); a_.grid(color="#e4e3df")
    fig.tight_layout(); fig.savefig(out / "axial_profiles.png", dpi=160); plt.close(fig)
    return summary


if __name__ == "__main__":
    main()
