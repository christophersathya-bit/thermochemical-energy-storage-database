"""Comparative figures for the TCES database -> figures/*.png

Run:  python scripts/plots.py
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

# Material group -> (colour, marker). Fixed order, never cycled.
GROUPS = {
    "Salt hydrates":  ("#2a78d6", "o"),
    "Hydroxides":     ("#eb6834", "s"),
    "Carbonates":     ("#1baf7a", "^"),
    "Metal oxides":   ("#eda100", "D"),
    "Metal hydrides": ("#e87ba4", "v"),
    "Ammoniates":     ("#008300", "P"),
    "Gas-phase":      ("#6250d6", "X"),
    "Sorption":       ("#e34948", "h"),
}
MECH2GROUP = {
    "Salt hydrate dehydration": "Salt hydrates", "Hydroxide dehydration": "Hydroxides",
    "Carbonate decarbonation": "Carbonates", "Metal oxide redox": "Metal oxides",
    "Metal hydride": "Metal hydrides", "Ammoniation": "Ammoniates",
    "Gas-phase reaction": "Gas-phase", "Physical adsorption": "Sorption", "Liquid absorption": "Sorption",
}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 10,
    "axes.titleweight": "bold", "axes.titlesize": 12, "axes.titlecolor": INK, "axes.titlelocation": "left",
})


def load():
    df = pd.read_csv(ROOT / "data" / "tces_database.csv")
    df["group"] = df["mechanism"].map(MECH2GROUP)
    df["short"] = df["reaction"].str.split(">").str[0].str.replace(r":\d+", "", regex=True).str.replace(r"\(.*?\)$", "", regex=True).str.strip()
    df.loc[df["family"] == "Sorption", "short"] = df["name"]
    # salt hydrates: "CuSO4 5→3" style labels so different hydration steps of one salt stay distinct
    import re as _re
    def _hyd(rx):
        left, right = rx.split(">")
        m1 = _re.match(r"\s*([A-Za-z0-9()]+)\.(\d*\.?\d*)H2O", left)
        if not m1:
            return None
        m2 = _re.match(r"\s*([A-Za-z0-9()]+)(?:\.(\d*\.?\d*)H2O)?", right)
        lo = (m2.group(2) or "1") if m2 and m2.group(2) is not None else "0"
        if m2 and m2.group(2) == "":
            lo = "1"
        return f"{m1.group(1)} {m1.group(2) or '1'}→{lo}"
    hyd = df["mechanism"] == "Salt hydrate dehydration"
    df.loc[hyd, "short"] = df.loc[hyd, "reaction"].map(lambda x: _hyd(x) or x)
    df["verified"] = ~df["verification"].fillna("").str.startswith("NOT verified")
    return df


def legend_handles(groups):
    return [Line2D([], [], marker=GROUPS[g][1], color="none", markerfacecolor=GROUPS[g][0],
                   markeredgecolor=SURF, markersize=9, label=g) for g in groups]


def save(fig, name):
    fig.savefig(FIG / name, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("saved", name)


def fig_temperature_vs_density(df):
    from adjustText import adjust_text
    d = df.dropna(subset=["ED_vol_GJ_m3"])
    fig, (ax, az) = plt.subplots(1, 2, figsize=(14, 6.2), gridspec_kw={"width_ratios": [1.6, 1]})
    for b in (200, 500):
        ax.axvline(b, color=INK2, lw=1, ls="--", alpha=.6)
    for x, t in [(100, "LOW"), (350, "MEDIUM"), (900, "HIGH")]:
        ax.text(x, 7.75, t, ha="center", color=INK2, fontsize=9, fontweight="bold")
    for a in (ax, az):
        for g, sub in d.groupby("group"):
            c, m = GROUPS[g]
            v, nv = sub[sub.verified], sub[~sub.verified]
            a.scatter(v["T_operating_C"], v["ED_vol_GJ_m3"], s=90, c=c, marker=m, edgecolors=SURF, linewidths=1.5, zorder=3)
            a.scatter(nv["T_operating_C"], nv["ED_vol_GJ_m3"], s=90, facecolors=SURF, marker=m, edgecolors=c, linewidths=1.8, zorder=3)
    # main panel: label only points outside the zoom box
    zoom = (d.T_operating_C < 200)
    texts = [ax.text(r.T_operating_C, r.ED_vol_GJ_m3, r.short, fontsize=8, color=INK2) for r in d[~zoom].itertuples()]
    adjust_text(texts, ax=ax, expand=(1.3, 1.6))
    ax.add_patch(plt.Rectangle((20, 0), 180, 3.1, fill=False, ec=INK2, lw=1, ls=":"))
    ax.text(130, 3.25, "zoom ->", color=INK2, fontsize=8, ha="center")
    ax.set_xlabel("Operating (charging) temperature  [°C]")
    ax.set_ylabel("Volumetric energy density  [GJ/m³ of material]")
    ax.set_title("Temperature vs. energy density of TCES systems")
    ax.set_ylim(0, 8); ax.set_xlim(0, 1400)
    ax.legend(handles=legend_handles([g for g in GROUPS if g in set(d.group)]), loc="upper left",
              bbox_to_anchor=(0.36, 0.98), frameon=False, ncol=2)
    texts = [az.text(r.T_operating_C, r.ED_vol_GJ_m3, r.short, fontsize=8, color=INK2) for r in d[zoom].itertuples()]
    az.set_xlim(20, 200); az.set_ylim(0, 3.1)
    az.text(0.99, 0.995, "hollow marker = value not yet verified against a source", transform=az.transAxes,
            ha="right", va="top", fontsize=7.5, color=INK2)
    adjust_text(texts, ax=az, expand=(1.4, 1.8), arrowprops=dict(arrowstyle="-", color=GRID, lw=0.8))
    az.set_title("Zoom: low-temperature systems (<200 °C)")
    az.set_xlabel("Operating (charging) temperature  [°C]")
    fig.tight_layout()
    save(fig, "01_temperature_vs_volumetric_density.png")


def fig_gravimetric_bar(df):
    d = df.dropna(subset=["ED_grav_MJ_kg"]).sort_values("ED_grav_MJ_kg")
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.barh(d["name"], d["ED_grav_MJ_kg"], color=[GROUPS[g][0] for g in d["group"]],
            edgecolor=SURF, linewidth=2, height=0.8)
    for i, v in enumerate(d["ED_grav_MJ_kg"]):
        ax.text(v + 0.05, i, f"{v:.2f}", va="center", fontsize=8, color=INK2)
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Gravimetric energy density  [MJ/kg]  (1 MJ/kg = 0.278 kWh/kg)")
    ax.set_title("Gravimetric energy density of all systems")
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUPS[g][0]) for g in GROUPS if g in set(d.group)]
    ax.legend(handles, [g for g in GROUPS if g in set(d.group)], loc="lower right", frameon=False)
    save(fig, "02_gravimetric_energy_density.png")


def fig_group_ranges(df):
    """Temperature window and energy density range of each material group (small multiples)."""
    order = [g for g in GROUPS if g in set(df.group)]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharey=True)
    for i, g in enumerate(order):
        sub = df[df.group == g]
        c = GROUPS[g][0]
        for ax, col in zip(axes, ["T_operating_C", "ED_grav_MJ_kg"]):
            v = sub[col].dropna()
            if v.empty:
                continue
            ax.plot([v.min(), v.max()], [i, i], color=c, lw=6, solid_capstyle="round", alpha=.35)
            ax.scatter(v, [i] * len(v), color=c, marker=GROUPS[g][1], s=60, edgecolors=SURF, linewidths=1.2, zorder=3)
    axes[0].set_yticks(range(len(order)), order)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Operating temperature  [°C]")
    axes[1].set_xlabel("Gravimetric energy density  [MJ/kg]")
    axes[0].set_title("Temperature window by material group")
    axes[1].set_title("Energy density range by material group")
    for ax in axes:
        ax.grid(axis="y", visible=False)
    save(fig, "03_group_ranges.png")


def fig_vant_hoff(df):
    """Equilibrium pressure lines ln(p/p0) = -dH/RT + dS/R for calc gas-solid reactions (per mol gas)."""
    R = 8.314
    d = df[(df.data_mode == "calc") & (df.subfamily == "Gas-solid")]
    groups = [g for g in GROUPS if g in set(d.group)]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), sharex=True, sharey=True)
    for ax, g in zip(axes.flat, groups):
        sub = d[d.group == g]
        for k, (_, r) in enumerate(sub.iterrows()):
            dh = r["dH_kJ_per_mol_gas"] * 1e3
            ds = r["dS_J_molK_rxn"] / r["nu_gas"]
            T = np.linspace(300, 1700, 300)
            p = np.exp(-dh / (R * T) + ds / R)
            ok = (p > 1e-4) & (p < 1e2)
            ax.semilogy(1000 / T[ok], p[ok], color=GROUPS[g][0], lw=2,
                        ls=["-", "--", ":", "-."][k % 4], label=f"{r['short']} ({r['T_eq_C']:.0f} °C @1 bar)")
        ax.legend(frameon=False, fontsize=7.5, loc="lower left" if g == "Salt hydrates" else "upper right")
        ax.axhline(1, color=INK2, lw=0.8, ls=":")
        ax.set_title(g, fontsize=10)
    for ax in axes.flat[len(groups):]:
        ax.set_visible(False)
    for ax in axes[-1]:
        ax.set_xlabel("1000 / T  [1/K]")
    for ax in axes[:, 0]:
        ax.set_ylabel("Equilibrium gas pressure [bar]")
    fig.suptitle("Van 't Hoff equilibrium lines (from standard data, per mol of gas)", x=0.01, ha="left", fontweight="bold")
    fig.tight_layout()
    save(fig, "04_vant_hoff_lines.png")


def fig_class_counts(df):
    t = df.groupby(["temperature_class", "group"]).size().unstack(fill_value=0)
    t = t.reindex(["Low (<200 C)", "Medium (200-500 C)", "High (>500 C)"])
    order = [g for g in GROUPS if g in t.columns]
    fig, ax = plt.subplots(figsize=(9, 4.5))
    left = np.zeros(len(t))
    for g in order:
        ax.barh(t.index, t[g], left=left, color=GROUPS[g][0], edgecolor=SURF, linewidth=2, label=g, height=0.6)
        left += t[g].values
    ax.invert_yaxis()
    ax.set_xlabel("Number of systems in the database")
    ax.set_title("Database coverage: systems per temperature class")
    ax.grid(axis="y", visible=False)
    ax.legend(frameon=False, bbox_to_anchor=(1, 1), loc="upper left")
    save(fig, "05_classification_counts.png")


if __name__ == "__main__":
    df = load()
    fig_temperature_vs_density(df)
    fig_gravimetric_bar(df)
    fig_group_ranges(df)
    fig_vant_hoff(df)
    fig_class_counts(df)
