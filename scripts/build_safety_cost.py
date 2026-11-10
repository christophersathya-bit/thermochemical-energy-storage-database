"""Join cost (Donkers 2017) and GHS hazard data (PubChem) to every database entry
-> data/safety_cost/material_safety_cost.csv (one row per id; empty = not sourced)."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "data" / "safety_cost"


def main():
    db = pd.read_csv(ROOT / "data" / "tces_database.csv")[["id", "name"]]
    mp = pd.read_csv(D / "material_compound_map.csv")
    ghs = pd.read_csv(D / "ghs_pubchem.csv").rename(columns={"compound": "hazard_compound"})
    cost = pd.read_csv(D / "cost_ld50_donkers2017.csv")
    g = mp.merge(ghs, on="hazard_compound", how="left")
    rank = {"Danger": 2, "Warning": 1}
    rows = []
    for i, sub in g.groupby("id"):
        st = sub[sub.role == "storage material"]
        wf = sub[sub.role == "working fluid"]
        s = st.iloc[0] if len(st) else None
        sig = max([x for x in sub.signal_word.dropna()], key=lambda x: rank.get(x, 0), default="")
        rows.append(dict(id=i,
                         hazard_compound=s.hazard_compound if s is not None else "",
                         pubchem_cid=int(s.pubchem_cid) if s is not None else "",
                         GHS_signal_word=sig,
                         GHS_pictograms=s.pictograms if s is not None and pd.notna(s.pictograms) else "",
                         GHS_H_statements=s.H_statements_majority if s is not None and pd.notna(s.H_statements_majority) else "",
                         working_fluid_hazard="; ".join(f"{r.hazard_compound}: {r.signal_word}, {r.pictograms}" for r in wf.itertuples())))
    h = pd.DataFrame(rows)
    out = db.merge(h, on="id", how="left").merge(cost[["id", "price_EUR_kg", "price_EUR_MJ", "LD50_mg_kg", "price_note"]], on="id", how="left")
    out["price_EUR_kWh"] = (out.price_EUR_MJ * 3.6).round(2)
    out["cost_source"] = out.price_EUR_kg.notna().map({True: "Donkers et al. 2017 Table 3 (2017 prices) [R13]", False: ""})
    out["hazard_source"] = out.hazard_compound.notna().map({True: "PubChem GHS Classification, retrieved 2026-11-11", False: ""})
    out["pubchem_cid"] = out.pubchem_cid.astype("Int64")
    out.to_csv(D / "material_safety_cost.csv", index=False)
    print(f"safety/cost: {out.GHS_signal_word.fillna('').ne('').sum()} entries with GHS signal word, "
          f"{out.price_EUR_kg.notna().sum()} with price, {out.LD50_mg_kg.notna().sum()} with LD50")


if __name__ == "__main__":
    main()
