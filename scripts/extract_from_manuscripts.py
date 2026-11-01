"""Semi-automatic data collection from manuscripts (PDF papers).

Put PDF papers in  manuscripts/  and run:
    python scripts/extract_from_manuscripts.py
Output:
    data/manuscript_candidates.csv  - one row per value found, with file, page and the
                                      sentence it came from, for MANUAL REVIEW.
    data/manuscript_tables/*.csv    - every table pdfplumber can detect.

What it looks for (regex, units normalised):
    reaction enthalpy      e.g. "ΔH = 104.4 kJ/mol", "enthalpy of 61 kJ mol-1"
    energy density         e.g. "2.8 GJ/m3", "628 kWh/m3", "1.4 MJ/kg", "1400 kJ/kg"
    temperatures           e.g. "dehydration at 120 °C", "T_eq = 505 °C"
    chemical formulas      salt hydrates / hydroxides / carbonates / oxides / hydrides near the value

Workflow: review the candidate rows, copy confirmed values into data/reactions_input.csv
(data_mode = lit, ref_ids = new reference in data/references.csv), then rerun
scripts/compute_properties.py and database/build_database.py.
Never add values to the database without checking them against the paper.
"""
from pathlib import Path
import re
import sys
import pandas as pd

try:
    import pdfplumber
except ImportError:
    sys.exit("pip install pdfplumber")

ROOT = Path(__file__).resolve().parents[1]
MS = ROOT / "manuscripts"
OUT_TABLES = ROOT / "data" / "manuscript_tables"

NUM = r"(\d+(?:[.,]\d+)?)"
PATTERNS = {
    "enthalpy_kJ_mol":   (rf"{NUM}\s*kJ\s*(?:/|·|\s)?\s*mol(?:-1|−1|⁻¹)?", 1.0),
    "ED_GJ_m3":          (rf"{NUM}\s*GJ\s*(?:/|·)?\s*m(?:3|³|-3|−3)", 1.0),
    "kWh_m3_to_ED_GJ_m3": (rf"{NUM}\s*kWh\s*(?:/|·)?\s*m(?:3|³|-3|−3)", 0.0036),
    "ED_MJ_kg":          (rf"{NUM}\s*MJ\s*(?:/|·)?\s*kg", 1.0),
    "kJ_kg_to_ED_MJ_kg": (rf"{NUM}\s*kJ\s*(?:/|·)?\s*kg", 0.001),
    "temperature_C":     (rf"{NUM}\s*(?:°|º|o)\s*C\b", 1.0),
}
FORMULA = re.compile(r"\b(?:[A-Z][a-z]?\d*(?:\([A-Z][a-z]?\d*\))?\d*){1,4}(?:[·.]\s?\d*(?:\.\d)?H2O|\(OH\)2|CO3|O\d|H2)\b")


def sentences(text):
    return re.split(r"(?<=[.;])\s+(?=[A-Z(])", text.replace("\n", " "))


def scan_pdf(pdf_path):
    rows = []
    with pdfplumber.open(pdf_path) as pdf:
        for pno, page in enumerate(pdf.pages, start=1):
            text = page.extract_text() or ""
            for sent in sentences(text):
                formulas = sorted(set(FORMULA.findall(sent)))
                for prop, (pat, factor) in PATTERNS.items():
                    for m in re.finditer(pat, sent, flags=re.I):
                        val = float(m.group(1).replace(",", ".")) * factor
                        rows.append(dict(file=pdf_path.name, page=pno, property=prop.split("_to_")[-1],
                                         value=round(val, 4), raw=m.group(0), formulas_in_sentence="; ".join(formulas),
                                         sentence=sent.strip()[:400]))
            for k, tbl in enumerate(page.extract_tables() or []):
                if tbl and len(tbl) > 1:
                    OUT_TABLES.mkdir(parents=True, exist_ok=True)
                    pd.DataFrame(tbl[1:], columns=[str(c) for c in tbl[0]]).to_csv(
                        OUT_TABLES / f"{pdf_path.stem}_p{pno}_t{k + 1}.csv", index=False)
    return rows


def main():
    pdfs = sorted(MS.glob("*.pdf"))
    if not pdfs:
        print(f"No PDFs found in {MS}. Copy your manuscripts there and rerun.")
        return
    rows = []
    for p in pdfs:
        print("scanning", p.name)
        rows += scan_pdf(p)
    df = pd.DataFrame(rows)
    df["verified"] = False
    out = ROOT / "data" / "manuscript_candidates.csv"
    df.to_csv(out, index=False)
    print(f"{len(df)} candidate values from {len(pdfs)} manuscripts -> {out.relative_to(ROOT)}")
    if len(df):
        print(df.groupby("property").size().to_string())


if __name__ == "__main__":
    main()
