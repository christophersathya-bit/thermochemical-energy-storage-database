"""ML surrogate of the 1D packed-bed reactor model (CFD + ML coupling).

1. Samples reactor design/operating conditions with a Latin hypercube.
2. Runs cfd/python_1d_reactor/reactor_model.py for every sample (Ca(OH)2 discharge).
3. Trains Gaussian-process and random-forest surrogates that predict
      - time to 90 % conversion  [min]
      - peak thermal power        [W]
      - maximum outlet temperature [°C]
   and reports test-set R2 + parity plots.
4. The trained surrogate evaluates a new design in milliseconds instead of seconds,
   which makes design optimisation / sensitivity studies cheap.

Usage
    python ml/reactor_surrogate.py --n 300          # generate + train (about 1-2 min)
    python ml/reactor_surrogate.py --reuse          # retrain on ml/results/surrogate_dataset.csv
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import qmc
from sklearn.ensemble import RandomForestRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel, ConstantKernel
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.model_selection import train_test_split
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "cfd" / "python_1d_reactor"))
from reactor_model import Params, run  # noqa: E402

OUT = ROOT / "ml" / "results"
OUT.mkdir(parents=True, exist_ok=True)

# design space: name -> (low, high)
SPACE = {
    "u_s": (0.05, 0.5),        # m/s
    "p_in": (0.2e5, 0.9e5),    # Pa  H2O partial pressure
    "T_in": (623.15, 703.15),  # K
    "h_w": (20.0, 300.0),      # W/m2/K
    "L": (0.1, 0.4),           # m
    "N0": (5000.0, 12000.0),   # mol/m3
}
TARGETS = {"t90_min": "Time to 90 % conversion [min]", "P_peak_W": "Peak thermal power [W]", "T_out_max_C": "Max outlet temperature [°C]"}


def generate(n, seed=0):
    sampler = qmc.LatinHypercube(d=len(SPACE), seed=seed)
    lo, hi = np.array([v[0] for v in SPACE.values()]), np.array([v[1] for v in SPACE.values()])
    X = qmc.scale(sampler.random(n), lo, hi)
    rows = []
    t0 = time.time()
    for i, x in enumerate(X):
        kw = dict(zip(SPACE, x))
        P = Params(name=f"lhs_{i}", T0=kw["T_in"], T_wall=kw["T_in"], t_end=3 * 3600, n_cells=50, n_out=241, **kw)
        try:
            k = run(P).kpis()
        except Exception as e:  # keep going if one sample fails
            print("sample", i, "failed:", e)
            continue
        rows.append({**kw, "t90_min": k["time_to_90pct_conversion_s"] / 60, "P_peak_W": k["peak_total_power_W"],
                     "T_out_max_C": k["T_out_max_C"]})
        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{n} simulations ({time.time() - t0:.0f} s)")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "surrogate_dataset.csv", index=False)
    return df


def train(df):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    df = df.dropna()
    X = df[list(SPACE)].values
    fig, axes = plt.subplots(1, len(TARGETS), figsize=(15, 4.6))
    report, models = [], {}
    for ax, (tgt, label) in zip(axes, TARGETS.items()):
        y = df[tgt].values
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.25, random_state=0)
        cands = {
            "Gaussian process": make_pipeline(StandardScaler(), GaussianProcessRegressor(
                ConstantKernel() * RBF(np.ones(X.shape[1])) + WhiteKernel(), normalize_y=True, n_restarts_optimizer=3, random_state=0)),
            "Random forest": RandomForestRegressor(n_estimators=500, random_state=0),
        }
        best = None
        for name, m in cands.items():
            m.fit(Xtr, ytr)
            yp = m.predict(Xte)
            r2, mae = r2_score(yte, yp), mean_absolute_error(yte, yp)
            report.append((tgt, name, r2, mae))
            if best is None or r2 > best[1]:
                best = (name, r2, m, yp)
        models[tgt] = best[2]
        ax.plot([y.min(), y.max()], [y.min(), y.max()], ls="--", color="#52514e", lw=1)
        ax.scatter(yte, best[3], s=40, color="#2a78d6", edgecolors="white")
        ax.set(xlabel=f"1D model: {label}", ylabel="Surrogate prediction", title=f"{best[0]}  R² = {best[1]:.3f}")
        ax.grid(color="#e4e3df"); ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("ML surrogate of the packed-bed reactor (test set)", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout(); fig.savefig(OUT / "surrogate_parity.png", dpi=160); plt.close(fig)
    rep = pd.DataFrame(report, columns=["target", "model", "R2_test", "MAE_test"])
    rep.to_csv(OUT / "surrogate_scores.csv", index=False)
    joblib.dump(dict(models=models, inputs=list(SPACE)), OUT / "surrogate_models.joblib")
    print(rep.round(3).to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--reuse", action="store_true")
    a = ap.parse_args()
    data = pd.read_csv(OUT / "surrogate_dataset.csv") if a.reuse else generate(a.n)
    train(data)
