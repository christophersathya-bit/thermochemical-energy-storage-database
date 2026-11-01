"""Build database/tces.db (SQLite) from the CSV files in data/.

Usage:
    python scripts/compute_properties.py      # refresh data/tces_database.csv
    python database/build_database.py         # build database/tces.db
    python database/build_database.py --query "SELECT name, T_operating_C FROM reactions WHERE ED_vol_GJ_m3 > 3"
"""
import argparse
import sqlite3
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DB = ROOT / "database" / "tces.db"


def build():
    DB.unlink(missing_ok=True)
    con = sqlite3.connect(DB)
    con.executescript((ROOT / "database" / "schema.sql").read_text())

    refs = pd.read_csv(ROOT / "data" / "references.csv")
    sp = pd.read_csv(ROOT / "data" / "thermo_species.csv")
    db = pd.read_csv(ROOT / "data" / "tces_database.csv")

    refs.to_sql("references_", con, if_exists="append", index=False)
    sp.to_sql("species", con, if_exists="append", index=False)
    db[["mechanism", "family", "subfamily"]].drop_duplicates().to_sql("classification", con, if_exists="append", index=False)
    cols = [c for c in db.columns if c not in ("family", "subfamily", "ref_ids")]
    db[cols].to_sql("reactions", con, if_exists="append", index=False)
    links = [(r.id, ref) for r in db.itertuples() for ref in str(r.ref_ids).split(";") if ref]
    con.executemany("INSERT INTO reaction_references VALUES (?, ?)", links)
    con.commit()
    n = con.execute("SELECT COUNT(*) FROM reactions").fetchone()[0]
    print(f"Built {DB.name}: {n} reactions, {len(sp)} species, {len(refs)} references")
    return con


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--query", help="SQL query to run after building")
    a = ap.parse_args()
    con = build()
    if a.query:
        print(pd.read_sql(a.query, con).to_string(index=False))
