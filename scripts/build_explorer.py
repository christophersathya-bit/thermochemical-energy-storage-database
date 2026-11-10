"""Build the interactive explorer (explorer/index.html) from the material tables.
Open explorer/index.html in a browser; it works offline except for the web fonts."""
from pathlib import Path
import datetime, json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEEP = ["id", "name", "family", "subfamily", "mechanism", "working_gas", "reaction", "density_discharged_g_cm3",
        "Cp298_discharged_J_kgK", "dH_kJ_mol_rxn", "dS_J_molK_rxn", "T_eq_C", "T_operating_C", "T_basis", "T_release_C",
        "ED_grav_MJ_kg", "ED_vol_GJ_m3", "ED_vol_basis", "ED_bed_GJ_m3", "volume_decrease_pct", "volume_change_calc_pct",
        "DRH_pct", "cyclability", "main_issue", "ref_ids", "verification", "source_location", "hazard_compound",
        "GHS_signal_word", "GHS_pictograms", "GHS_H_statements", "working_fluid_hazard", "price_EUR_kg", "price_EUR_kWh",
        "LD50_mg_kg"]


def main():
    m = pd.read_csv(ROOT / "data" / "thermophysical" / "materials_properties.csv")[KEEP]
    notes = pd.read_csv(ROOT / "data" / "tces_database.csv")[["id", "notes"]]
    m = m.merge(notes, on="id", how="left").replace({np.nan: None})
    rows = [{k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if v is not None} for r in m.to_dict("records")]
    r = pd.read_csv(ROOT / "data" / "references.csv").replace({np.nan: None})
    refs = {x.ref_id: dict(cite=f"{x.authors} ({x.year}) {x.title}. {x.journal or ''} {x.volume_pages or ''}".strip(),
                           url=(f"https://doi.org/{x.doi}" if x.doi else (x.url or ""))) for x in r.itertuples()}
    t = (ROOT / "explorer" / "template.html").read_text(encoding="utf-8")
    body = (t.replace("__DATA__", json.dumps(rows, ensure_ascii=False, separators=(",", ":")))
             .replace("__REFS__", json.dumps(refs, ensure_ascii=False, separators=(",", ":")))
             .replace("__DATE__", "2026-11-11"))
    (ROOT / "explorer" / "explorer_body.html").write_text(body, encoding="utf-8")
    full = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
            '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
            + body.replace("<style>", "<style>\nbody{margin:0}\n", 1).split("<div class=\"wrap\">")[0]
            + '</head>\n<body>\n<div class="wrap">' + body.split('<div class="wrap">', 1)[1] + "\n</body>\n</html>\n")
    (ROOT / "explorer" / "index.html").write_text(full, encoding="utf-8")
    print(f"explorer: {len(rows)} materials, {len(refs)} references -> explorer/index.html")


if __name__ == "__main__":
    main()
