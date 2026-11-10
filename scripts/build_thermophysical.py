"""Write the material and thermophysical property CSV files -> data/thermophysical/

  materials_properties.csv    one row per storage material (all 65): thermodynamic, energy-storage and
                              thermophysical properties (molar mass, densities, heat capacities, conductivity)
  species_properties.csv      one row per chemical species: M, density, formation enthalpy, entropy, Cp(298 K)
  species_cp_vs_T.csv         Cp, S, H-H298 vs temperature from NIST-JANAF (copied from the tables)
  gas_properties.csv          humid-N2 / steam properties vs T and vapour fraction (cfd/properties.py)
  class_property_ranges.csv   property ranges per material class (Kiyabu et al. 2026)
  cfd_bed_properties.csv      the bed / solid properties used by both CFD models

Every value carries a source column. Empty cells mean "no sourced value yet" - nothing is filled in by guess.
"""
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "thermophysical"
sys.path.insert(0, str(ROOT / "cfd"))
from properties import gas_mixture, BedProps  # noqa: E402

# Densities of highest / lowest hydrate printed in Donkers et al. 2017, Table 3 (p.52) [R13]
DONKERS_RHO = {  # id: (rho_highest_hydrate, rho_lowest_hydrate)
    "S06": (2.18, 2.33), "S07": (1.56, 1.89), "S08": (1.58, 1.86), "S10": (1.76, 2.09), "S11": (2.67, 3.51),
    "S12": (1.84, 2.15), "S13": (2.23, 2.58), "S14": (2.2, 2.78), "S15": (2.51, 3.39),
}
K_SALT_HYDRATES = "0.3-1.3 W/m/K for salt hydrates (class range, Kiyabu et al. 2026 p.2103)"


def parse(rx):
    if ">" not in rx or ":" not in rx:
        return None, None
    left, right = rx.split(">")
    side = lambda s: [(t.split(":")[0].strip(), float(t.split(":")[1])) for t in s.split("+")]
    return side(left), side(right)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    db = pd.read_csv(ROOT / "data" / "tces_database.csv")
    sp = pd.read_csv(ROOT / "data" / "thermo_species.csv").set_index("species")
    cpT = pd.read_csv(OUT / "species_cp_vs_T.csv")
    cp298 = cpT[cpT.T_K == 298.15].set_index("species")

    # ---- species_properties.csv ------------------------------------------------
    s = sp.reset_index()[["species", "phase", "molar_mass_g_mol", "density_g_cm3", "density_source",
                          "dHf_kJ_mol", "S_J_molK", "source", "source_url"]].copy()
    s["Cp298_J_molK"] = s.species.map(cp298["Cp_J_molK"])
    s["Cp298_J_kgK"] = s.species.map(cp298["Cp_J_kgK"])
    s["Cp_source"] = np.where(s.Cp298_J_molK.notna(), "NIST-JANAF (same table as dHf, S)", "")
    s = s.rename(columns={"source": "dHf_S_source", "source_url": "dHf_S_Cp_url"})
    s.to_csv(OUT / "species_properties.csv", index=False)

    # ---- materials_properties.csv ----------------------------------------------
    rows = []
    for _, r in db.iterrows():
        left, right = parse(r.reaction)
        d = dict(id=r.id, name=r["name"], family=r.family, subfamily=r.subfamily, mechanism=r.mechanism,
                 working_gas=r.working_gas, reaction=r.reaction)
        dis = chg = None
        if left:
            sl = [x for x, _ in left if sp.loc[x, "phase"] == "s"]
            sr = [x for x, _ in right if sp.loc[x, "phase"] == "s"]
            dis = sl[0] if sl else None
            chg = sr[0] if sr else None
        d["discharged_species"] = dis or ""
        d["charged_species"] = chg or ""
        d["M_discharged_g_mol"] = sp.loc[dis, "molar_mass_g_mol"] if dis else r.M_left_g_mol
        d["M_charged_g_mol"] = sp.loc[chg, "molar_mass_g_mol"] if chg else np.nan
        if r.id in DONKERS_RHO:
            d["density_discharged_g_cm3"], d["density_charged_g_cm3"] = DONKERS_RHO[r.id]
            d["density_source"] = "Donkers et al. 2017 Table 3 (highest / lowest hydrate) [R13]"
        else:
            d["density_discharged_g_cm3"] = sp.loc[dis, "density_g_cm3"] if dis else r.density_left_g_cm3
            d["density_charged_g_cm3"] = sp.loc[chg, "density_g_cm3"] if chg else np.nan
            d["density_source"] = ("CRC handbook crystal density (not yet re-fetched)"
                                   if (dis or pd.notna(d["density_discharged_g_cm3"])) else "")
        for tag, x in (("discharged", dis), ("charged", chg)):
            d[f"Cp298_{tag}_J_kgK"] = cp298.loc[x, "Cp_J_kgK"] if x in cp298.index else np.nan
        d["Cp_source"] = "NIST-JANAF 298.15 K" if any(x in cp298.index for x in (dis, chg) if x) else ""
        d["thermal_conductivity_W_mK"] = np.nan
        d["thermal_conductivity_note"] = K_SALT_HYDRATES if r.mechanism == "Salt hydrate dehydration" else ""
        for c in ["dH_kJ_mol_rxn", "dS_J_molK_rxn", "nu_gas", "dH_kJ_per_mol_gas", "T_eq_C", "T_operating_C", "T_basis",
                  "T_release_C", "temperature_class", "ED_grav_MJ_kg", "ED_grav_kWh_kg", "ED_vol_GJ_m3", "ED_vol_basis",
                  "ED_bed_GJ_m3", "volume_decrease_pct", "DRH_pct", "cyclability", "main_issue", "ref_ids",
                  "verification", "source_location"]:
            d[c] = r[c]
        if pd.notna(d["density_discharged_g_cm3"]) and pd.notna(d["density_charged_g_cm3"]) and pd.isna(d["volume_decrease_pct"]) and left:
            # crystal volume change from molar masses and densities (only where both densities exist)
            v_dis = d["M_discharged_g_mol"] / d["density_discharged_g_cm3"]
            n_chg = [n for x, n in right if x == chg][0] / [n for x, n in left if x == dis][0]
            v_chg = n_chg * d["M_charged_g_mol"] / d["density_charged_g_cm3"]
            d["volume_change_calc_pct"] = round(100 * (v_chg - v_dis) / v_dis, 1)
        else:
            d["volume_change_calc_pct"] = np.nan
        rows.append(d)
    m = pd.DataFrame(rows)
    sc = ROOT / "data" / "safety_cost" / "material_safety_cost.csv"
    if sc.exists():  # cost (Donkers 2017) and GHS hazards (PubChem), see data/safety_cost/README.md
        m = m.merge(pd.read_csv(sc).drop(columns=["name"]), on="id", how="left")
    m.to_csv(OUT / "materials_properties.csv", index=False)

    # ---- gas_properties.csv ----------------------------------------------------
    g = []
    for T in (293.15, 373.15, 473.15, 573.15, 673.15, 773.15, 873.15):
        for y in (0.0, 0.01, 0.1, 0.25, 0.5, 1.0):
            if T < 373.15 and y > 0.02:
                continue  # steam partial pressure above saturation at 1 bar
            p = gas_mixture(T, 1e5, y * 1e5)
            g.append(dict(T_C=round(T - 273.15), y_H2O=y, p_bar=1.0, M_g_mol=round(p.M * 1e3, 3), rho_kg_m3=round(p.rho, 4),
                          cp_J_kgK=round(p.cp, 1), mu_Pa_s=f"{p.mu:.4e}", k_W_mK=round(p.k, 5), Pr=round(p.Pr, 3),
                          method="NIST Shomate cp; Sutherland mu, k; Wilke mixing (cfd/properties.py)"))
    pd.DataFrame(g).to_csv(OUT / "gas_properties.csv", index=False)

    # ---- class_property_ranges.csv ---------------------------------------------
    c = pd.read_csv(ROOT / "data" / "kiyabu2026" / "table4_class_overview.csv")
    c["thermal_conductivity"] = np.where(c.material == "Salt hydrates", "0.3-1.3 W/m/K (p.2103)", "")
    c["source"] = "Kiyabu et al. 2026, Mater. Horiz. 13, 2087 (Table 4, p.2097)"
    c.to_csv(OUT / "class_property_ranges.csv", index=False)

    # ---- cfd_bed_properties.csv ------------------------------------------------
    b = BedProps()
    gp = gas_mixture(673.15, 1e5, 0.5e5)
    d, f = b.ergun(gp.mu)
    rows = [("bed porosity", b.eps, "-", "assumed design value (same in 1D and OpenFOAM)"),
            ("granule diameter", b.dp, "m", "assumed design value"),
            ("granule apparent density", b.rho_particle, "kg/m3", "assumed (porous Ca(OH)2 granules, internal porosity ~0.5)"),
            ("solid heat capacity", b.cp_solid, "J/kg/K", "assumed average; JANAF Ca(OH)2 1180 J/kg/K at 298 K, 1328 at 400 K"),
            ("reactive solid per bed volume N0", round(b.N0, 1), "mol/m3", "(1-eps) rho_p / M"),
            ("bed volumetric heat capacity", round(b.rho_cp_bed, 0), "J/m3/K", "(1-eps) rho_p cp_s"),
            ("effective bed conductivity", 0.4, "W/m/K", "assumed; within 0.3-1.3 W/m/K salt-hydrate range (Kiyabu 2026) - fit to experiments"),
            ("Darcy coefficient d", f"{d:.4e}", "1/m2", "Ergun: 150(1-eps)^2/(eps^3 dp^2)"),
            ("Forchheimer coefficient f", f"{f:.4e}", "1/m", "Ergun: 3.5(1-eps)/(eps^3 dp)"),
            ("reaction enthalpy", 106.8, "kJ/mol H2O", "Schaube et al. 2012 equilibrium fit [R5]"),
            ("reaction entropy", 137.2, "J/mol/K", "Schaube et al. 2012 equilibrium fit [R5]")]
    pd.DataFrame(rows, columns=["property", "value", "unit", "source_or_basis"]).to_csv(OUT / "cfd_bed_properties.csv", index=False)

    print(f"materials {len(m)} | species {len(s)} | Cp(T) rows {len(cpT)} | gas rows {len(g)}")
    print("filled per column (materials):")
    print(m.notna().sum().loc[["density_discharged_g_cm3", "density_charged_g_cm3", "Cp298_discharged_J_kgK",
                               "Cp298_charged_J_kgK", "ED_grav_MJ_kg", "ED_vol_GJ_m3", "dH_kJ_mol_rxn", "T_operating_C"]].to_string())


if __name__ == "__main__":
    main()
