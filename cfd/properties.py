"""Material properties shared by the 1D reactor model and the OpenFOAM case.

Gas phase: ideal mixture of a reactive gas (H2O) and an inert carrier (N2).
  cp   : NIST Chemistry WebBook Shomate equations (per species), mole-fraction mixing
  mu   : Sutherland law per species, Wilke mixing rule
  k    : Sutherland law per species, Wilke-type (Mason-Saxena) mixing
  rho  : ideal gas
Solid / bed: Ca(OH)2 granules (defaults), all derived from ONE porosity so both CFD
models use the same bed.

Run `python properties.py` to print the property table used by the default case.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

R = 8.314462618

# --- species data ---------------------------------------------------------------
M = {"N2": 28.0134e-3, "H2O": 18.0153e-3, "CO2": 44.0095e-3, "O2": 31.9988e-3}  # kg/mol

# NIST Shomate coefficients  cp = A + B t + C t^2 + D t^3 + E / t^2   [J/mol/K], t = T/1000
SHOMATE = {
    "N2": [((100, 500), (28.98641, 1.853978, -9.647459, 16.63537, 0.000117)),
           ((500, 2000), (19.50583, 19.88705, -8.598535, 1.369784, 0.527601))],
    "H2O": [((500, 1700), (30.09200, 6.832514, 6.793435, -2.534480, 0.082139))],
}
# Sutherland: mu = mu0 (T/T0)^1.5 (T0+S)/(T+S);  same form for k
SUTH_MU = {"N2": (1.663e-5, 273.15, 107.0), "H2O": (1.12e-5, 350.0, 1064.0)}
SUTH_K = {"N2": (0.0242, 273.15, 150.0), "H2O": (0.0181, 300.0, 2200.0)}


def cp_molar(sp: str, T: float) -> float:
    for (lo, hi), (A, B, C, D, E) in SHOMATE[sp]:
        if lo <= T <= hi:
            t = T / 1000
            return A + B * t + C * t ** 2 + D * t ** 3 + E / t ** 2
    if sp == "H2O" and T < 500:          # ideal-gas steam below the Shomate range
        return 33.58 + (cp_molar("H2O", 500) - 33.58) * (T - 298.15) / (500 - 298.15)
    raise ValueError(f"T={T} K outside cp range for {sp}")


def _suth(c, T):
    v0, T0, S = c
    return v0 * (T / T0) ** 1.5 * (T0 + S) / (T + S)


def _wilke(x, prop, T):
    sps = list(x)
    out = 0.0
    for i in sps:
        den = 0.0
        for j in sps:
            phi = (1 + (prop[i] / prop[j]) ** 0.5 * (M[j] / M[i]) ** 0.25) ** 2 / (8 * (1 + M[i] / M[j])) ** 0.5
            den += x[j] * phi
        out += x[i] * prop[i] / den
    return out


@dataclass
class GasProps:
    T: float
    p: float
    y_H2O: float
    M: float        # kg/mol
    rho: float      # kg/m3
    cp: float       # J/kg/K
    mu: float       # Pa s
    k: float        # W/m/K

    @property
    def Pr(self):
        return self.cp * self.mu / self.k

    @property
    def rho_cp(self):
        return self.rho * self.cp


def gas_mixture(T: float, p_tot: float = 1e5, p_H2O: float = 0.5e5) -> GasProps:
    y = p_H2O / p_tot
    x = {"H2O": y, "N2": 1 - y}
    x = {s: v for s, v in x.items() if v > 0}
    Mmix = sum(x[s] * M[s] for s in x)
    cp = sum(x[s] * cp_molar(s, T) for s in x) / Mmix
    mu = _wilke(x, {s: _suth(SUTH_MU[s], T) for s in x}, T)
    k = _wilke(x, {s: _suth(SUTH_K[s], T) for s in x}, T)
    return GasProps(T, p_tot, y, Mmix, p_tot * Mmix / (R * T), cp, mu, k)


@dataclass
class BedProps:
    """Packed bed of porous Ca(OH)2 granules (defaults)."""
    eps: float = 0.40              # inter-particle bed porosity (used by BOTH models)
    dp: float = 2.0e-3             # granule diameter [m]
    rho_particle: float = 1100.0   # apparent granule density [kg/m3] (granule internal porosity ~0.5)
    M_solid: float = 74.093e-3     # Ca(OH)2 [kg/mol]
    cp_solid: float = 1300.0       # J/kg/K  (Ca(OH)2 / CaO average at 400-500 C)

    @property
    def N0(self):                  # mol of reactive solid per m3 of bed
        return (1 - self.eps) * self.rho_particle / self.M_solid

    @property
    def rho_cp_bed(self):          # J/m3/K
        return (1 - self.eps) * self.rho_particle * self.cp_solid

    def ergun(self, mu):
        d = 150 * (1 - self.eps) ** 2 / (self.eps ** 3 * self.dp ** 2)
        f = 3.5 * (1 - self.eps) / (self.eps ** 3 * self.dp)
        return d, f

    def dp_per_length(self, u, gas: GasProps):
        d, f = self.ergun(gas.mu)
        return gas.mu * d * u + 0.5 * gas.rho * f * u * u


if __name__ == "__main__":
    print(f"{'T [C]':>6} {'y_H2O':>6} {'M g/mol':>8} {'rho':>7} {'cp':>7} {'mu e5':>7} {'k':>7} {'Pr':>5}")
    for T in (293.15, 673.15, 823.15):
        for y in (0.01, 0.5):
            g = gas_mixture(T, 1e5, y * 1e5)
            print(f"{T - 273.15:6.0f} {y:6.2f} {g.M * 1e3:8.2f} {g.rho:7.3f} {g.cp:7.0f} {g.mu * 1e5:7.3f} {g.k:7.4f} {g.Pr:5.2f}")
    b = BedProps()
    print(f"\nBed: eps={b.eps}, dp={b.dp * 1e3} mm, N0={b.N0:.0f} mol/m3, rho*cp={b.rho_cp_bed:.3e} J/m3/K")
