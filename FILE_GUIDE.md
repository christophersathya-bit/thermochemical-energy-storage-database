# File guide: what every file in this repository is

## ⭐ The material database (start here)

| File | What it is |
|---|---|
| [`data/tces_database.csv`](data/tces_database.csv) | **The thermochemical material database.** 80 materials, one row each: reaction, classification, ΔH, ΔS, temperatures, energy densities, cyclability, references, verification status |
| [`data/thermophysical/materials_properties.csv`](data/thermophysical/materials_properties.csv) | **The same 80 materials + thermophysical properties** (molar mass, densities, heat capacity, conductivity note). Use this if you want everything in one file |
| [`database/tces.db`](database/tces.db) | The same data as an SQLite database (open with DB Browser for SQLite, Python, etc.) |

Open a CSV file on GitHub to see it as a table; use **Download raw file** (top right of the file view) to save it and open in Excel.

## `data/` – inputs, sources and property tables

| File | What it is |
|---|---|
| `data/reactions_input.csv` | Input sheet the database is built from (edit this to add materials) |
| `data/thermo_species.csv` | ΔfH°, S°, molar mass, density of 38 chemical species, with NIST-JANAF links |
| `data/references.csv` | Full citation + DOI for every reference code (R1–R20) used in `ref_ids` |
| `data/classification_summary.csv` | Number of materials, temperature and energy-density range per mechanism |
| `data/verification_report.csv` | Cross-checks of database values against NIST and papers |
| `data/thermophysical/species_properties.csv` | 38 species: M, density, ΔfH°, S°, Cp at 298 K |
| `data/thermophysical/species_cp_vs_T.csv` | Cp, S, H−H298 at 298–1200 K (NIST-JANAF) |
| `data/thermophysical/gas_properties.csv` | N₂ / steam density, cp, viscosity, conductivity, Pr at 20–600 °C |
| `data/thermophysical/class_property_ranges.csv` | Property ranges per material class |
| `data/thermophysical/cfd_bed_properties.csv` | Bed / solid properties used in the CFD |
| `data/thermophysical/README.md` | Column definitions and data gaps of the property files |
| `data/ntsoukpoe2014/`, `andre2016/`, `yan2015/`, `schaube2012/`, `block2016/`, `sogutoglu2018/`, `literature_extra/` | Tables copied from ScienceDirect papers (one folder per paper, each with a README giving the DOI) |
| `data/kiyabu2026/*.csv` | Tables 1–5, 7, 8 copied from Kiyabu et al., Mater. Horiz. 2026 (targets, prototypes, cooling pairs, class overview, salts, liquid pairs, composites) |

## `database/` – SQLite database

| File | What it is |
|---|---|
| `database/tces.db` | SQLite database (tables `reactions`, `species`, `references_`, `classification`, `kiyabu2026_*`; view `v_reactions_full`) |
| `database/schema.sql` | Table definitions |
| `database/build_database.py` | Builds `tces.db` from the CSV files; `--query "SQL"` runs a query |

## `figures/` – plots used in the README

| File | What it shows |
|---|---|
| `00_classification_tree.png` / `.mmd` | How TCES materials are classified |
| `01_temperature_vs_volumetric_density.png` | Temperature vs energy density of all materials |
| `02_gravimetric_energy_density.png` | Ranking by energy per kg |
| `03_group_ranges.png` | Temperature and energy range per material group |
| `04_vant_hoff_lines.png` | Equilibrium pressure–temperature lines |
| `05_classification_counts.png` | Number of materials per temperature class |
| `06_1d_reactor_caoh2_discharge.png` | 1D reactor simulation result |
| `07_openfoam_bed_fields.png`, `08_openfoam_axial_profiles.png` | OpenFOAM temperature / pressure results |
| `09`–`11_ml_*.png` | Machine-learning results |
| `12_cfd_validation.png` | CFD validation checks |
| `13_composites_and_liquid_pairs.png` | Composites and liquid absorption pairs |

## `docs/`

| File | What it is |
|---|---|
| `docs/DATA_VERIFICATION.md` | How every entry was verified, with cross-checks |

## `cfd/` – CFD models

| File | What it is |
|---|---|
| `cfd/properties.py` | Gas and bed material properties shared by both CFD models |
| `cfd/validate_cfd.py` | Runs the 19 validation checks |
| `cfd/VALIDATION.md`, `cfd/validation.json` | Validation results (19/19 passed) |
| `cfd/python_1d_reactor/reactor_model.py` | 1D transient packed-bed reactor model (Python) |
| `cfd/python_1d_reactor/cases/*.json` | Input cases: Ca(OH)₂ charge / discharge, MgSO₄ discharge |
| `cfd/python_1d_reactor/results/*/` | Results per case: plot, `fields.csv`, `outlet.csv`, `kpis.json` |
| `cfd/openfoam/packedBedReactor/` | OpenFOAM case (`0/`, `constant/`, `system/`, `Allrun`, `Allclean`) |
| `cfd/openfoam/setup_from_1d.py` | Writes properties and conditions into the OpenFOAM case |
| `cfd/openfoam/postprocess.py` | Plots and checks the OpenFOAM results |

## `ml/` – machine learning

| File | What it is |
|---|---|
| `ml/tces_ml.py` | Turning-temperature and energy-density models, material classifier, screening |
| `ml/reactor_surrogate.py` | ML surrogate trained on 300 reactor simulations |
| `ml/results/*` | Model scores, parity plots, screening ranking, surrogate dataset |

## `scripts/` – processing scripts

| File | What it does |
|---|---|
| `scripts/compute_properties.py` | Calculates ΔH, ΔS, T_eq, energy densities → `data/tces_database.csv` |
| `scripts/build_thermophysical.py` | Writes the files in `data/thermophysical/` |
| `scripts/classification.py` | Classification tree and summary |
| `scripts/plots.py` | Comparative figures |
| `scripts/verification_report.py` | Verification report |
| `scripts/extract_from_manuscripts.py` | Extracts candidate values from PDF papers in `manuscripts/` |

## Other files

| File | What it is |
|---|---|
| `README.md` | Overview, figures, results |
| `run_all.py` | Rebuilds everything (`python run_all.py`) |
| `requirements.txt` | Python packages needed |
| `manuscripts/README.md` | Put PDF papers here for extraction |
| `.github/workflows/ci.yml` | Automatic rebuild check on GitHub |
| `LICENSE` | MIT (code), CC BY 4.0 (data) |
