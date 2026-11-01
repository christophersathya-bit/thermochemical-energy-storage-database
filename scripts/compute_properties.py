"""
Compute thermodynamic properties of every TCES reaction and write the master table
data/tces_database.csv.

Two data modes (column `data_mode` in data/reactions_input.csv):
  calc : reaction enthalpy/entropy computed from standard formation data
         (data/thermo_species.csv, 298.15 K) -> dH, dS, T_eq = dH/dS (1 bar gas),
         gravimetric and volumetric energy density.
  lit  : values taken directly from published manuscripts (refs in references.csv).

Energy densities are given per kg / per m3 of the *discharged* (energy-poor) solid,
i.e. the left-hand side of the charging reaction  A(s) -> B(s) + gas.
For gas-phase reactions the basis is the total mass of the feed.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def parse_side(side: str):
    out = []
    for term in side.split("+"):
        sp, nu = term.strip().split(":")
        out.append((sp.strip(), float(nu)))
    return out


def temperature_class(t):
    if pd.isna(t):
        return np.nan
    if t < 200:
        return "Low (<200 C)"
    if t <= 500:
        return "Medium (200-500 C)"
    return "High (>500 C)"


def main():
    sp = pd.read_csv(DATA / "thermo_species.csv").set_index("species")
    rx = pd.read_csv(DATA / "reactions_input.csv")
    rows = []
    for _, r in rx.iterrows():
        d = r.to_dict()
        if r["data_mode"] == "calc":
            left, right = (parse_side(s) for s in r["reaction"].split(">"))
            dH = sum(n * sp.loc[s, "dHf_kJ_mol"] for s, n in right) - sum(n * sp.loc[s, "dHf_kJ_mol"] for s, n in left)
            dS = sum(n * sp.loc[s, "S_J_molK"] for s, n in right) - sum(n * sp.loc[s, "S_J_molK"] for s, n in left)
            gas_nu = sum(n for s, n in right if sp.loc[s, "phase"] == "g") - sum(n for s, n in left if sp.loc[s, "phase"] == "g")
            solids_left = [(s, n) for s, n in left if sp.loc[s, "phase"] == "s"]
            basis = solids_left if solids_left else left
            mass = sum(n * sp.loc[s, "molar_mass_g_mol"] for s, n in basis)
            rho = sp.loc[solids_left[0][0], "density_g_cm3"] if solids_left else np.nan
            T_eq_K = dH * 1000 / dS
            nbs = [sname for sname, _ in left + right if not str(sp.loc[sname, "source"]).startswith("NIST-JANAF")]
            d["verification"] = ("computed from NIST-JANAF" if not nbs else
                                 "computed; NBS 1982 values not yet re-fetched for " + ", ".join(sorted(set(nbs))))
            d.update(
                dH_kJ_mol_rxn=dH, dS_J_molK_rxn=dS, nu_gas=gas_nu,
                dH_kJ_per_mol_gas=dH / gas_nu if gas_nu else np.nan,
                T_eq_C=T_eq_K - 273.15,
                T_operating_C=r["T_charge_C_lit"] if pd.notna(r["T_charge_C_lit"]) else T_eq_K - 273.15,
                M_left_g_mol=mass, density_left_g_cm3=rho,
                ED_grav_MJ_kg=dH / mass,
                ED_vol_GJ_m3=dH / mass * rho if not np.isnan(rho) else np.nan,
                ED_grav_kWh_kg=dH / mass / 3.6,
            )
        else:  # literature values
            dH = r["dH_kJ_mol_lit"]
            M = r["M_left_g_mol_lit"]
            rho = r["density_left_g_cm3_lit"]
            edg = r["ED_grav_MJ_kg_lit"] if pd.notna(r["ED_grav_MJ_kg_lit"]) else (dH / M if pd.notna(dH) and pd.notna(M) else np.nan)
            edv = r["ED_vol_GJ_m3_lit"] if pd.notna(r["ED_vol_GJ_m3_lit"]) else (edg * rho if pd.notna(rho) else np.nan)
            d.update(
                dH_kJ_mol_rxn=dH, dS_J_molK_rxn=r["dS_J_molK_lit"], nu_gas=np.nan,
                dH_kJ_per_mol_gas=np.nan, T_eq_C=np.nan, T_operating_C=r["T_charge_C_lit"],
                M_left_g_mol=M, density_left_g_cm3=rho,
                ED_grav_MJ_kg=edg, ED_vol_GJ_m3=edv, ED_grav_kWh_kg=edg / 3.6 if pd.notna(edg) else np.nan,
            )
        if r["data_mode"] == "calc":
            d["T_basis"] = "T_eq at 1 bar reactive gas (298 K data)" if pd.isna(r["T_charge_C_lit"]) else "typical operating T"
        elif "20 mbar" in str(r["notes"]):
            d["T_basis"] = "dehydration T at 20 mbar water vapour"
        elif "storage (charging)" in str(r["notes"]):
            d["T_basis"] = "storage (charging) T of prototype / cycle"
        elif "heat release" in str(r["notes"]) or "T_operating = T_release" in str(r["notes"]):
            d["T_basis"] = "heat-release T (charging T not reported)"
        else:
            d["T_basis"] = "typical charging T from source"
        d["temperature_class"] = temperature_class(d["T_operating_C"])
        rows.append(d)

    out = pd.DataFrame(rows)
    keep = ["id", "name", "reaction", "family", "subfamily", "mechanism", "working_gas", "data_mode",
            "dH_kJ_mol_rxn", "dS_J_molK_rxn", "nu_gas", "dH_kJ_per_mol_gas", "T_eq_C", "T_operating_C", "T_basis",
            "temperature_class", "M_left_g_mol", "density_left_g_cm3", "ED_grav_MJ_kg", "ED_grav_kWh_kg",
            "ED_vol_GJ_m3", "ED_vol_basis", "ED_bed_GJ_m3", "T_release_C", "volume_decrease_pct", "DRH_pct", "cyclability", "main_issue", "ref_ids", "verification", "source_location", "notes"]
    out = out[keep].round(3)
    out.to_csv(DATA / "tces_database.csv", index=False)
    print(out[["id", "name", "T_operating_C", "ED_grav_MJ_kg", "ED_vol_GJ_m3", "temperature_class"]].to_string())


if __name__ == "__main__":
    main()
