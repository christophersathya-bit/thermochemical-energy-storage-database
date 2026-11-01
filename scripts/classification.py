"""Build the TCES classification tree from the database.

Writes:
  figures/00_classification_tree.mmd   (Mermaid source; GitHub renders it in Markdown)
  data/classification_summary.csv      (counts and ranges per branch)
Render to PNG (optional):  mmdc -i figures/00_classification_tree.mmd -o figures/00_classification_tree.png -w 2000

Classification logic used in this database
  Level 1  Family      : Sorption (physical bond)  |  Chemical reaction (chemical bond)
  Level 2  Subfamily   : Adsorption / Absorption   |  Gas-solid / Gas-gas (homogeneous)
  Level 3  Mechanism   : hydration, carbonation, redox, hydride, ammoniation, ...
  Level 4  Material    : individual reaction pair
  Cross-cut            : temperature class  Low <200 C | Medium 200-500 C | High >500 C
"""
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def node_id(s):
    return re.sub(r"[^A-Za-z0-9]", "_", s)


def main():
    df = pd.read_csv(ROOT / "data" / "tces_database.csv")
    lines = ["flowchart LR", '  TCES["Thermochemical<br/>energy storage"]']
    style = {"Sorption": "#fde2e2", "Chemical reaction": "#dbe8fa"}
    for fam, f in df.groupby("family", sort=False):
        fid = node_id(fam)
        lines.append(f'  TCES --> {fid}["{fam}"]')
        for sub, s in f.groupby("subfamily", sort=False):
            sid = node_id(fam + sub)
            lines.append(f'  {fid} --> {sid}["{sub}"]')
            for mech, m in s.groupby("mechanism", sort=False):
                mid = node_id(mech)
                t = f"{m.T_operating_C.min():.0f}-{m.T_operating_C.max():.0f} °C" if len(m) > 1 else f"{m.T_operating_C.iloc[0]:.0f} °C"
                n = f"{len(m)} system" + ("s" if len(m) > 1 else "")
                lines.append(f'  {sid} --> {mid}["{mech}<br/><small>{n} · {t}</small>"]')
                lines.append(f"  style {mid} fill:{style[fam]},stroke:#999")
    (ROOT / "figures").mkdir(exist_ok=True)
    (ROOT / "figures" / "00_classification_tree.mmd").write_text("\n".join(lines) + "\n")

    summ = (df.groupby(["family", "subfamily", "mechanism"])
              .agg(n_systems=("id", "count"),
                   T_min_C=("T_operating_C", "min"), T_max_C=("T_operating_C", "max"),
                   ED_grav_min_MJ_kg=("ED_grav_MJ_kg", "min"), ED_grav_max_MJ_kg=("ED_grav_MJ_kg", "max"),
                   ED_vol_max_GJ_m3=("ED_vol_GJ_m3", "max"))
              .round(2).reset_index())
    summ.to_csv(ROOT / "data" / "classification_summary.csv", index=False)
    print(summ.to_string(index=False))


if __name__ == "__main__":
    main()
