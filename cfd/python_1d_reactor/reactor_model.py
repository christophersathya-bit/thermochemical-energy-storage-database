"""
1D transient packed-bed reactor model for gas-solid thermochemical energy storage.

Finite-volume (method of lines) model of a fixed bed through which a reactive gas
(H2O, CO2, O2, ...) diluted in an inert carrier flows.  Default case: Ca(OH)2/CaO
discharge (hydration) at the conditions of a lab-scale reactor.

Governing equations (x = axial coordinate, local thermal equilibrium gas/solid):

  Reactive-gas balance (mole fraction y, constant total pressure, inert carrier conserved)
      eps * c_tot(T) * dy/dt + d/dx[ F_I * y/(1-y) ] = eps * D_ax * c_tot * d2y/dx2 - nu * N0 * dX/dt
      F_I = inert molar flux (constant), c_tot = p_tot/(R T), partial pressure p = y * p_tot
      (heating the gas at constant total pressure does NOT raise the vapour pressure)
  Energy balance (bed + gas)
      (rho*cp)_eff * dT/dt + rho_g*cp_g*u_s * dT/dx = lambda_eff * d2T/dx2
                              + N0 * dH * dX/dt - h_w * (4/D) * (T - T_wall)
  Solid conversion (X = fraction of the discharged / energy-poor form, e.g. Ca(OH)2)
      p_eq(T) = p0 * exp(-dH/(R*T) + dS/R)                  (van 't Hoff, per mol gas)
      discharge (p > p_eq):  dX/dt =  k_d(T) * (1 - X) * (p/p_eq - 1)
      charge    (p < p_eq):  dX/dt = -k_c(T) *  X      * (1 - p/p_eq)
      k(T) = A * exp(-Ea / (R*T))

Boundary conditions: Dirichlet inlet (y_in = p_in/p_tot, T_in), zero-gradient outlet.
Numerics: first-order upwind convection, central diffusion, stiff BDF integration.

Run:
    python reactor_model.py                       # default Ca(OH)2 discharge
    python reactor_model.py --config cases/caoh2_charge.json
    python reactor_model.py --config cases/mgso4_discharge.json --vtk
Outputs go to results/<case name>/ : CSV fields, outlet history, PNG plots, optional VTK.
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import diags, bmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from properties import gas_mixture  # noqa: E402

R = 8.314462618  # J/mol/K
_trapz = getattr(np, "trapezoid", None) or np.trapz  # numpy 1.x / 2.x


@dataclass
class Params:
    name: str = "caoh2_discharge"
    mode: str = "discharge"            # 'discharge' (exothermic) or 'charge' (endothermic)
    # reaction (per mol of gas)
    dH: float = 106.8e3                # J/mol gas  Ca(OH)2 -> CaO + H2O, Schaube 2012 equilibrium fit
    dS: float = 137.2                  # J/mol/K    ln(p/bar) = -12845/T + 16.508  ->  T_eq(1 bar) = 505 C
    nu: float = 1.0                    # mol gas per mol reaction
    A_d: float = 1.3945e4              # 1/s  discharge pre-exponential   (order of magnitude, verify)
    Ea_d: float = 89.5e3               # J/mol
    A_c: float = 1.9425e12             # 1/s  charge pre-exponential      (order of magnitude, verify)
    Ea_c: float = 191.9e3              # J/mol
    # bed
    L: float = 0.20                    # m bed length
    D: float = 0.05                    # m bed diameter
    eps: float = 0.40                  # inter-particle bed porosity (same as OpenFOAM case)
    dp: float = 2.0e-3                 # m granule diameter (Ergun pressure drop)
    N0: float = 8908.0                 # mol/m3 bed = (1-eps)*rho_particle/M = 0.6*1100/0.0741  (cfd/properties.py)
    X0: float = 0.0                    # initial conversion (0 = fully charged)
    rho_cp_bed: float = 8.58e5         # J/m3/K  (1-eps)*rho_particle*cp_s = 0.6*1100*1300
    lam_eff: float = 0.4               # W/m/K   effective axial conductivity
    D_ax: float = 1.0e-4               # m2/s    axial dispersion
    h_w: float = 150.0                 # W/m2/K  wall / heat-exchanger coefficient (0 = adiabatic)
    T_wall: float = 673.15             # K       heat-transfer-fluid temperature (indirect heat exchange)
    # inlet / initial state
    u_s: float = 0.20                  # m/s superficial velocity
    p_tot: float = 1.0e5               # Pa
    p_in: float = 0.5e5                # Pa partial pressure of reactive gas at inlet
    T_in: float = 673.15               # K inlet temperature
    T0: float = 673.15                 # K initial bed temperature
    p0_bed: float = 0.0                # Pa initial reactive gas partial pressure in bed
    rho_cp_gas: Optional[float] = None # J/m3/K rho_g*cp_g; None = humid-N2 mixture from cfd/properties.py at inlet
    # numerics
    n_cells: int = 80
    t_end: float = 3600.0              # s
    n_out: int = 121

    @classmethod
    def from_json(cls, path):
        return cls(**json.loads(Path(path).read_text()))


class PackedBed:
    def __init__(self, p: Params):
        self.p = p
        self.gas = gas_mixture(p.T_in, p.p_tot, p.p_in)
        if p.rho_cp_gas is None:
            p.rho_cp_gas = self.gas.rho_cp
        n = p.n_cells
        self.dx = p.L / n
        self.x = (np.arange(n) + 0.5) * self.dx

    def p_eq(self, T):
        return self.p.p_tot * np.exp(-self.p.dH / (R * T) + self.p.dS / R)

    def rate(self, X, yg, T):
        """dX/dt (1/s); positive = discharge (gas uptake, heat release)."""
        P = self.p
        pg = np.clip(yg, 0.0, 1.0) * P.p_tot
        peq = self.p_eq(T)
        ratio = pg / peq
        rd = P.A_d * np.exp(-P.Ea_d / (R * T)) * (1 - X) * np.clip(ratio - 1, 0, None)
        rc = P.A_c * np.exp(-P.Ea_c / (R * T)) * X * np.clip(1 - ratio, 0, None)
        return rd - rc

    def _transport(self, f, f_in, u, diff):
        dx = self.dx
        fl = np.concatenate(([f_in], f[:-1]))
        conv = -u * (f - fl) / dx
        fpad = np.concatenate(([2 * f_in - f[0]], f, [f[-1]]))  # ghost cells
        dif = diff * (fpad[2:] - 2 * f + fpad[:-2]) / dx ** 2
        return conv + dif

    def rhs(self, t, y):
        P, n = self.p, self.p.n_cells
        X, yg, T = y[:n], y[n:2 * n], y[2 * n:]
        r = self.rate(X, yg, T)
        y_in = P.p_in / P.p_tot
        ctot = P.p_tot / (R * T)
        F_I = P.u_s * P.p_tot / (R * P.T_in) * (1 - y_in)          # inert molar flux, mol/m2/s
        yc = np.clip(yg, 0.0, 0.999)
        F = np.concatenate(([F_I * y_in / (1 - y_in)], F_I * yc / (1 - yc)))   # upwind face fluxes
        ypad = np.concatenate(([2 * y_in - yg[0]], yg, [yg[-1]]))
        disp = P.eps * P.D_ax * ctot * (ypad[2:] - 2 * yg + ypad[:-2]) / self.dx ** 2
        dc = (-(F[1:] - F[:-1]) / self.dx + disp - P.nu * P.N0 * r) / (P.eps * ctot)
        dT = (self._transport(T, P.T_in, P.rho_cp_gas * P.u_s, P.lam_eff)
              + P.N0 * P.dH * r * P.nu - P.h_w * 4 / P.D * (T - P.T_wall)) / P.rho_cp_bed
        return np.concatenate((r, dc, dT))

    def jac_sparsity(self):
        n = self.p.n_cells
        tri = diags([1.0, 1.0, 1.0], [-1, 0, 1], shape=(n, n))
        I = diags([1.0], [0], shape=(n, n))
        return bmat([[I, I, I], [I, tri, I], [I, I, tri]]).tocsr()

    def solve(self):
        P, n = self.p, self.p.n_cells
        y0 = np.concatenate((np.full(n, P.X0), np.full(n, P.p0_bed / P.p_tot), np.full(n, P.T0)))
        t_eval = np.linspace(0, P.t_end, P.n_out)
        sol = solve_ivp(self.rhs, (0, P.t_end), y0, method="BDF", t_eval=t_eval,
                        jac_sparsity=self.jac_sparsity(), rtol=1e-5, atol=1e-7)
        if not sol.success:
            raise RuntimeError(sol.message)
        X, c, T = sol.y[:n], sol.y[n:2 * n], sol.y[2 * n:]
        return Result(self, sol.t, X, c, T)


@dataclass
class Result:
    model: PackedBed
    t: np.ndarray
    X: np.ndarray
    c: np.ndarray        # mole fraction of reactive gas
    T: np.ndarray

    def outlet_power_W(self):
        """Thermal power carried by the gas: (rho cp)_g * u_s * A * (T_out - T_in)."""
        P = self.model.p
        A = np.pi * P.D ** 2 / 4
        return P.rho_cp_gas * P.u_s * A * (self.T[-1] - P.T_in)

    def wall_power_W(self):
        """Heat transferred to the heat-transfer fluid through the wall (positive = heat delivered)."""
        P = self.model.p
        A = np.pi * P.D ** 2 / 4
        return (P.h_w * 4 / P.D * (self.T - P.T_wall)).sum(axis=0) * A * self.model.dx

    def kpis(self):
        P = self.model.p
        A = np.pi * P.D ** 2 / 4
        Xm = self.X.mean(axis=0)
        stored = P.N0 * P.dH * P.nu * A * P.L  # J at full conversion
        target = 0.9 if P.mode == "discharge" else 0.1
        reached = np.where((Xm >= target) if P.mode == "discharge" else (Xm <= target))[0]
        if len(reached) and reached[0] > 0:            # linear interpolation of the crossing time
            j = reached[0]
            t90 = float(np.interp(target, sorted([Xm[j - 1], Xm[j]]),
                                  [self.t[j - 1], self.t[j]] if Xm[j] >= Xm[j - 1] else [self.t[j], self.t[j - 1]]))
        else:
            t90 = np.nan
        pw = self.outlet_power_W() + self.wall_power_W()
        # global energy balance:  Q_reaction = Q_gas + Q_wall + Q_inlet_conduction + dU_sensible
        dx = self.model.dx
        q_rxn = P.N0 * P.dH * P.nu * A * dx * (self.X[:, -1] - self.X[:, 0]).sum()
        p_cond_in = P.lam_eff * A * (self.T[0] - P.T_in) / (dx / 2)
        q_out = _trapz(pw + p_cond_in, self.t)
        dU = P.rho_cp_bed * A * dx * (self.T[:, -1] - self.T[:, 0]).sum()
        g = self.model.gas
        dp_ergun = (g.mu * 150 * (1 - P.eps) ** 2 / (P.eps ** 3 * P.dp ** 2) * P.u_s
                    + 1.75 * g.rho * (1 - P.eps) / (P.eps ** 3 * P.dp) * P.u_s ** 2) * P.L
        return dict(
            energy_balance_error_pct=float(100 * (q_rxn - q_out - dU) / abs(q_rxn)) if q_rxn else 0.0,
            pressure_drop_ergun_Pa=float(dp_ergun),
            gas_rho_kg_m3=g.rho, gas_cp_J_kgK=g.cp, gas_mu_Pa_s=g.mu, gas_k_W_mK=g.k,
            final_mean_conversion=float(Xm[-1]),
            time_to_90pct_conversion_s=t90,
            T_out_max_C=float(self.T[-1].max() - 273.15),
            T_out_min_C=float(self.T[-1].min() - 273.15),
            peak_total_power_W=float(np.abs(pw).max()),
            mean_total_power_W=float(np.abs(_trapz(pw, self.t)) / self.t[-1]),
            mean_gas_power_W=float(np.abs(_trapz(self.outlet_power_W(), self.t)) / self.t[-1]),
            reaction_energy_kJ=float(abs(Xm[-1] - P.X0) * stored / 1e3),
            capacity_kJ=float(stored / 1e3),
        )

    def save(self, outdir: Path, vtk=False):
        import pandas as pd
        outdir.mkdir(parents=True, exist_ok=True)
        P, x = self.model.p, self.model.x
        rows = []
        for k, tk in enumerate(self.t):
            for i, xi in enumerate(x):
                rows.append((tk, xi, self.X[i, k], self.c[i, k] * P.p_tot, self.T[i, k] - 273.15))
        pd.DataFrame(rows, columns=["t_s", "x_m", "X", "p_gas_Pa", "T_C"]).to_csv(outdir / "fields.csv", index=False)
        pd.DataFrame(dict(t_s=self.t, T_out_C=self.T[-1] - 273.15, X_mean=self.X.mean(axis=0),
                          gas_power_W=self.outlet_power_W(), wall_power_W=self.wall_power_W())).to_csv(outdir / "outlet.csv", index=False)
        (outdir / "kpis.json").write_text(json.dumps(dict(params=asdict(P), kpis=self.kpis()), indent=2))
        if vtk:
            self._write_vtk(outdir / "vtk")
        self.plot(outdir)

    def _write_vtk(self, d: Path):
        """Legacy VTK rectilinear grid per time step (open the series in ParaView)."""
        d.mkdir(exist_ok=True)
        x = self.model.x
        for k, tk in enumerate(self.t):
            with open(d / f"bed_{k:04d}.vtk", "w") as f:
                f.write(f"# vtk DataFile Version 3.0\nTCES bed t={tk:.1f}s\nASCII\nDATASET RECTILINEAR_GRID\n")
                f.write(f"DIMENSIONS {len(x)} 1 1\nX_COORDINATES {len(x)} float\n{' '.join(f'{v:.5f}' for v in x)}\n")
                f.write("Y_COORDINATES 1 float\n0\nZ_COORDINATES 1 float\n0\n")
                f.write(f"POINT_DATA {len(x)}\n")
                for name, arr in (("T_C", self.T[:, k] - 273.15), ("X", self.X[:, k]), ("p_gas_Pa", self.c[:, k] * self.model.p.p_tot)):
                    f.write(f"SCALARS {name} float 1\nLOOKUP_TABLE default\n{' '.join(f'{v:.5g}' for v in arr)}\n")

    def plot(self, outdir: Path):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        P, x = self.model.p, self.model.x * 100
        fig, ax = plt.subplots(2, 2, figsize=(11, 7.5))
        idx = np.linspace(0, len(self.t) - 1, 6).astype(int)
        cmap = plt.cm.Blues(np.linspace(0.35, 1, len(idx)))
        for col, k in zip(cmap, idx):
            ax[0, 0].plot(x, self.T[:, k] - 273.15, color=col, lw=2, label=f"{self.t[k] / 60:.0f} min")
            ax[0, 1].plot(x, self.X[:, k], color=col, lw=2)
        ax[0, 0].set(xlabel="Bed position [cm]", ylabel="Temperature [°C]", title="Temperature profiles")
        ax[0, 0].legend(frameon=False, fontsize=8, loc="best", ncol=2)
        ax[0, 1].set(xlabel="Bed position [cm]", ylabel="Conversion X [-]", title="Reaction front (X = discharged fraction)", ylim=(-0.02, 1.02))
        ax[1, 0].plot(self.t / 60, self.T[-1] - 273.15, color="#2a78d6", lw=2, label="outlet")
        ax[1, 0].axhline(P.T_in - 273.15, color="#52514e", ls=":", lw=1, label="inlet")
        Teq = P.dH / (P.dS - R * np.log(P.p_in / P.p_tot)) - 273.15
        ax[1, 0].axhline(Teq, color="#eb6834", ls="--", lw=1, label=f"T_eq at p_in = {Teq:.0f} °C")
        ax[1, 0].set(xlabel="Time [min]", ylabel="Temperature [°C]", title="Outlet temperature")
        ax[1, 0].legend(frameon=False, fontsize=8)
        ax[1, 1].plot(self.t / 60, self.wall_power_W(), color="#2a78d6", lw=2, label="to heat-transfer fluid (wall)")
        ax[1, 1].plot(self.t / 60, self.outlet_power_W(), color="#eb6834", lw=2, ls="--", label="carried by gas (out - in)")
        ax[1, 1].legend(frameon=False, fontsize=8)
        ax[1, 1].set(xlabel="Time [min]", ylabel="Thermal power [W]", title="Thermal power")
        for a in ax.flat:
            a.grid(color="#e4e3df"); a.spines[["top", "right"]].set_visible(False)
        fig.suptitle(f"1D packed-bed TCES reactor - {P.name}", x=0.01, ha="left", fontweight="bold")
        fig.tight_layout()
        fig.savefig(outdir / f"{P.name}.png", dpi=160)
        plt.close(fig)


def run(params: Params, outdir: Path | None = None, vtk=False):
    res = PackedBed(params).solve()
    if outdir is not None:
        res.save(outdir, vtk=vtk)
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", help="JSON file with Params fields")
    ap.add_argument("--vtk", action="store_true", help="also write VTK files for ParaView")
    a = ap.parse_args()
    P = Params.from_json(a.config) if a.config else Params()
    here = Path(__file__).resolve().parent
    res = run(P, here / "results" / P.name, vtk=a.vtk)
    print(json.dumps(res.kpis(), indent=2))
