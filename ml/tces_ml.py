"""Machine learning on the TCES database.

Tasks
  1. turning-temperature : predict T_eq (1 bar) from dH per mol gas, gas type, molar mass, density
                           (the model has to learn dS, which is roughly constant per working gas)
  2. energy-density      : predict gravimetric energy density from simple descriptors
  3. classify            : predict the material group (salt hydrate, carbonate, ...) from
                           operating temperature, energy density and working gas
  4. screen              : rank systems for a target application temperature (multi-criteria score)

The database is small (tens of systems), so every model is evaluated with leave-one-out
cross-validation and the numbers should be read as indicative. Add data with
scripts/extract_from_manuscripts.py to make the models stronger.

Usage
    python ml/tces_ml.py all
    python ml/tces_ml.py screen --T 150 --tol 60
"""
import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, r2_score, accuracy_score, confusion_matrix
from sklearn.model_selection import LeaveOneOut, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plots import MECH2GROUP, GROUPS, INK2, GRID  # noqa: E402  (shared colours / groups)

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", module="sklearn")
OUT = ROOT / "ml" / "results"
OUT.mkdir(parents=True, exist_ok=True)
GASES = ["H2O", "CO2", "O2 (air)", "H2", "NH3"]


def load():
    df = pd.read_csv(ROOT / "data" / "tces_database.csv")
    df["group"] = df["mechanism"].map(MECH2GROUP)
    for g in GASES:
        df[f"gas_{g}"] = (df["working_gas"] == g).astype(int)
    df["log_M"] = np.log(df["M_left_g_mol"])
    return df


def _plt():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def parity(y, yhat, groups, title, unit, fname):
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.5, 5.2))
    lo, hi = min(y.min(), yhat.min()), max(y.max(), yhat.max())
    pad = 0.05 * (hi - lo)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK2, lw=1, ls="--")
    for g in GROUPS:
        m = groups == g
        if m.any():
            ax.scatter(y[m], yhat[m], s=60, c=GROUPS[g][0], marker=GROUPS[g][1], edgecolors="white", lw=1, label=g, zorder=3)
    ax.set(xlabel=f"Reference {unit}", ylabel=f"ML prediction (leave-one-out) {unit}", title=title)
    ax.grid(color=GRID); ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / fname, dpi=170); plt.close(fig)


def loo_regression(models, X, y):
    res = {}
    for name, mdl in models.items():
        yhat = cross_val_predict(mdl, X, y, cv=LeaveOneOut())
        res[name] = (yhat, mean_absolute_error(y, yhat), r2_score(y, yhat))
    return res


def task_turning_temperature(df):
    d = df[(df.data_mode == "calc") & (df.subfamily == "Gas-solid")].copy()
    feats = ["dH_kJ_per_mol_gas", "log_M", "density_left_g_cm3"] + [f"gas_{g}" for g in GASES]
    X, y = d[feats].values, d["T_eq_C"].values
    models = {
        "Linear regression": make_pipeline(StandardScaler(), LinearRegression()),
        "Random forest": RandomForestRegressor(n_estimators=400, random_state=0),
        "Gaussian process": make_pipeline(StandardScaler(), GaussianProcessRegressor(
            ConstantKernel() * RBF(length_scale=np.ones(len(feats))) + WhiteKernel(), normalize_y=True, random_state=0)),
    }
    res = loo_regression(models, X, y)
    best = min(res, key=lambda k: res[k][1])
    parity(y, res[best][0], d.group.values, f"Turning temperature - {best}", "[°C]", "parity_turning_temperature.png")
    return pd.DataFrame([(k, "T_eq [°C]", v[1], v[2]) for k, v in res.items()], columns=["model", "target", "MAE", "R2"])


def task_energy_density(df):
    d = df.dropna(subset=["dH_kJ_mol_rxn", "M_left_g_mol", "ED_grav_MJ_kg"]).copy()
    d["dH_log"] = np.log(d["dH_kJ_mol_rxn"].abs())
    feats = ["dH_log", "log_M", "T_operating_C"]
    X, y = d[feats].values, d["ED_grav_MJ_kg"].values
    models = {
        "Linear regression": make_pipeline(StandardScaler(), LinearRegression()),
        "Random forest": RandomForestRegressor(n_estimators=400, random_state=0),
    }
    res = loo_regression(models, X, y)
    best = min(res, key=lambda k: res[k][1])
    parity(y, res[best][0], d.group.values, f"Gravimetric energy density - {best}", "[MJ/kg]", "parity_energy_density.png")
    return pd.DataFrame([(k, "ED_grav [MJ/kg]", v[1], v[2]) for k, v in res.items()], columns=["model", "target", "MAE", "R2"])


def task_classify(df):
    d = df.dropna(subset=["T_operating_C", "ED_grav_MJ_kg"]).copy()
    feats = ["T_operating_C", "ED_grav_MJ_kg"] + [f"gas_{g}" for g in GASES]
    X, y = d[feats].values, d["group"].values
    clf = RandomForestClassifier(n_estimators=500, random_state=0, class_weight="balanced")
    yhat = cross_val_predict(clf, X, y, cv=LeaveOneOut())
    acc = accuracy_score(y, yhat)
    labels = [g for g in GROUPS if g in set(y)]
    cm = confusion_matrix(y, yhat, labels=labels)
    plt = _plt()
    fig, ax = plt.subplots(figsize=(7, 6))
    ax.imshow(cm, cmap="Blues")
    ax.set_xticks(range(len(labels)), labels, rotation=45, ha="right"); ax.set_yticks(range(len(labels)), labels)
    for i in range(len(labels)):
        for j in range(len(labels)):
            if cm[i, j]:
                ax.text(j, i, cm[i, j], ha="center", va="center", color="white" if cm[i, j] > cm.max() / 2 else "black")
    ax.set(xlabel="Predicted group", ylabel="True group", title=f"Material-group classifier (LOO accuracy {acc:.0%})")
    fig.tight_layout(); fig.savefig(OUT / "confusion_matrix_classifier.png", dpi=170); plt.close(fig)
    clf.fit(X, y)
    imp = pd.Series(clf.feature_importances_, index=feats).sort_values(ascending=False)
    imp.to_csv(OUT / "classifier_feature_importance.csv", header=["importance"])
    return pd.DataFrame([("Random forest", "material group (accuracy)", np.nan, acc)], columns=["model", "target", "MAE", "R2"])


CYC = {"very good": 1.0, "good": 0.8, "moderate": 0.5, "unknown": 0.3, "poor": 0.2}


def task_screen(df, T, tol, w_vol=0.4, w_grav=0.2, w_cyc=0.3, w_T=0.1):
    """Score = weighted sum of normalised criteria for systems that operate within T +/- tol."""
    d = df[(df.T_operating_C - T).abs() <= tol].copy()
    if d.empty:
        print(f"No system within {T} +/- {tol} °C"); return d
    norm = lambda s: (s - s.min()) / (s.max() - s.min()) if s.max() > s.min() else s * 0 + 1
    d["s_vol"] = norm(d["ED_vol_GJ_m3"].fillna(0))
    d["s_grav"] = norm(d["ED_grav_MJ_kg"].fillna(0))
    d["s_cyc"] = d["cyclability"].map(CYC).fillna(0.3)
    d["s_T"] = 1 - (d.T_operating_C - T).abs() / tol
    d["score"] = w_vol * d.s_vol + w_grav * d.s_grav + w_cyc * d.s_cyc + w_T * d.s_T
    d = d.sort_values("score", ascending=False)
    cols = ["id", "name", "group", "T_operating_C", "ED_vol_GJ_m3", "ED_grav_MJ_kg", "cyclability", "main_issue", "score"]
    d[cols].round(3).to_csv(OUT / f"screening_T{T:g}.csv", index=False)
    print(d[cols].round(2).to_string(index=False))
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("task", choices=["all", "turning-temperature", "energy-density", "classify", "screen"])
    ap.add_argument("--T", type=float, default=150, help="target temperature for screening [°C]")
    ap.add_argument("--tol", type=float, default=60, help="temperature window +/- [°C]")
    a = ap.parse_args()
    df = load()
    rows = []
    if a.task in ("all", "turning-temperature"):
        rows.append(task_turning_temperature(df))
    if a.task in ("all", "energy-density"):
        rows.append(task_energy_density(df))
    if a.task in ("all", "classify"):
        rows.append(task_classify(df))
    if rows:
        summary = pd.concat(rows, ignore_index=True).rename(columns={"R2": "R2 / accuracy"})
        summary.to_csv(OUT / "model_summary.csv", index=False)
        print(summary.round(3).to_string(index=False))
    if a.task in ("all", "screen"):
        print(f"\nScreening for {a.T} +/- {a.tol} °C")
        task_screen(df, a.T, a.tol)


if __name__ == "__main__":
    main()
