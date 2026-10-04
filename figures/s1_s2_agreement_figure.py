"""
S1 vs S2 plume agreement figure (three panels).

  (a) S1 plume area vs S2 plume area
  (b) distance between the S1 and S2 plume centres
  (c) overlap between the S1 and S2 plume masks (IoU)

Input:  matched_pairs_with_areas.csv, written by 07_measure_plumes.py
Output: fig_s1_s2_agreement.png, saved next to this script
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------
if Path("/home/jovyan/s2").exists():          # TIDE server
    CSV = Path("/home/jovyan/s2/07_figures/matched_pairs_with_areas.csv")
else:                                          # Mac external drive
    CSV = Path("/Volumes/External/TJ/010_optical/07_figures/matched_pairs_with_areas.csv")
OUT_PNG = Path(__file__).resolve().parent / "fig_s1_s2_agreement.png"

BLUE = "steelblue"         # both sensors detected a plume
ORANGE = "darkorange"      # only S1 detected a plume
GREEN = "mediumseagreen"   # only S2 detected a plume
GREY = "0.75"              # only one sensor detected a plume (panel c)
RED = "firebrick"          # best-fit line and median

# ---------------------------------------------------------------------------
# Load the pairs
# ---------------------------------------------------------------------------
pairs = pd.read_csv(CSV)
pairs = pairs[pairs["S2_Rating"] == 3]                              # cloud-free S2 only
pairs = pairs.dropna(subset=["S1_Area_km2", "S2_Area_km2"])         # mask from both sensors

s1_found = pairs["S1_Area_km2"] > 0
s2_found = pairs["S2_Area_km2"] > 0
both = pairs[s1_found & s2_found]           # both sensors detected a plume
s1_only = pairs[s1_found & ~s2_found]       # S2 mask is empty
s2_only = pairs[~s1_found & s2_found]       # S1 mask is empty
n_one = len(s1_only) + len(s2_only)

print(f"{len(pairs)} pairs: {len(both)} both, {len(s1_only)} S1 only, {len(s2_only)} S2 only")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def stats_box(ax, text, y, va):
    """White box with summary numbers in the right-hand corner of a panel."""
    ax.text(0.97, y, text, transform=ax.transAxes, ha="right", va=va, fontsize=9,
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))


def panel_title(ax, title, subtitle):
    """Bold title, with a smaller grey line underneath saying which pairs are shown."""
    ax.text(0, 1.10, title, transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.text(0, 1.03, subtitle, transform=ax.transAxes, fontsize=9.5, color="0.35")


def box_and_dots(ax, values, colors, xlabel, unit):
    """Horizontal box plot with every pair drawn as a dot on top."""
    values = np.asarray(values, dtype=float)
    ax.boxplot(values, vert=False, positions=[0], widths=0.5, whis=(0, 100),
               showfliers=False, patch_artist=True,
               boxprops=dict(facecolor="#dce6f0"), medianprops=dict(color=RED, lw=2))

    # spread the dots up and down a little so they don't sit on top of each other
    jitter = np.random.default_rng(0).uniform(-0.18, 0.18, len(values))
    ax.scatter(values, jitter, s=40, c=colors, edgecolor="k", lw=0.6, zorder=3, clip_on=False)

    q1, median, q3 = np.percentile(values, [25, 50, 75])
    stats_box(ax, f"n = {len(values)} pairs\n"
                  f"median = {median:.1f} {unit}\n"
                  f"IQR = {q1:.1f}–{q3:.1f} {unit}\n"
                  f"range = {values.min():.1f}–{values.max():.1f} {unit}", y=0.04, va="bottom")

    ax.set_xlim(left=0)
    ax.set_ylim(-1.9, 0.6)                  # box near the top, room for the stats below
    ax.set_yticks([])
    ax.set_xlabel(xlabel)
    ax.spines[["left", "right", "top"]].set_visible(False)
    ax.grid(axis="x", color="0.9")
    ax.set_box_aspect(1)                    # square, same size as panel (a)


fig, (ax_a, ax_b, ax_c) = plt.subplots(1, 3, figsize=(21, 7))
fig.subplots_adjust(wspace=0.32)

# ---------------------------------------------------------------------------
# (a) Plume area: S1 vs S2
# ---------------------------------------------------------------------------
for group, marker, color, label in [(both, "o", BLUE, "both sensors detected plume"),
                                    (s1_only, "v", ORANGE, "S1 only (S2 area = 0)"),
                                    (s2_only, "s", GREEN, "S2 only (S1 area = 0)")]:
    ax_a.scatter(group["S1_Area_km2"], group["S2_Area_km2"], marker=marker, s=45, color=color,
                 edgecolor="k", lw=0.6, zorder=3, clip_on=False,
                 label=f"{label} (n = {len(group)})")

x, y = pairs["S1_Area_km2"], pairs["S2_Area_km2"]
top = max(x.max(), y.max()) * 1.05
ax_a.plot([0, top], [0, top], "--", color="gray", lw=1, label="1:1 line")
ax_a.text(top * 0.50, top * 0.66, "S2 larger", color="0.45", rotation=45, ha="center")
ax_a.text(top * 0.66, top * 0.50, "S1 larger", color="0.45", rotation=45, ha="center")

slope, intercept = np.polyfit(x, y, 1)
r2 = np.corrcoef(x, y)[0, 1] ** 2
fit_x = np.linspace(0, top, 200)
fit_y = intercept + slope * fit_x
ax_a.plot(fit_x[fit_y >= 0], fit_y[fit_y >= 0], color=RED, lw=1.3, label="least-squares fit")

stats_box(ax_a, f"$R^2$ = {r2:.3f}  (n = {len(pairs)})\n"
                f"y = {slope:.3f}x + {intercept:.2f}\n"
                f"S1 larger in {(x > y).sum()} of {len(pairs)} pairs\n"
                f"S2 larger in {(y > x).sum()} of {len(pairs)} pairs", y=0.97, va="top")

ax_a.set_xlim(0, top)
ax_a.set_ylim(0, top)
ax_a.set_box_aspect(1)
ax_a.set_xlabel("S1 plume area (km²)")
ax_a.set_ylabel("S2 plume area (km²)")
ax_a.legend(fontsize=8, loc="upper left")
panel_title(ax_a, "(a) Plume area", f"All pairs (n = {len(pairs)})")

# ---------------------------------------------------------------------------
# (b) Centroid distance: only possible when both masks contain a plume
# ---------------------------------------------------------------------------
distance_km = both["Centroid_Dist_m"] / 1000
box_and_dots(ax_b, distance_km, BLUE, "Distance between S1 and S2 plume centroids (km)", "km")
panel_title(ax_b, "(b) Centroid distance",
            f"Pairs where both sensors detected a plume (n = {len(both)})")

# ---------------------------------------------------------------------------
# (c) Mask overlap (IoU): pairs where only one sensor saw a plume have IoU = 0
# ---------------------------------------------------------------------------
one_sensor = ~(s1_found & s2_found)
colors = np.where(one_sensor, GREY, BLUE)
box_and_dots(ax_c, pairs["IoU_percent"], colors,
             "Mask overlap, intersection over union (%)", "%")

ax_c.scatter([], [], s=40, color=BLUE, edgecolor="k",
             label=f"both sensors detected plume (n = {len(both)})")
ax_c.scatter([], [], s=40, color=GREY, edgecolor="k",
             label=f"one sensor only, IoU = 0 (n = {n_one})")
ax_c.legend(fontsize=8, loc="center left", bbox_to_anchor=(0, 0.45))
panel_title(ax_c, "(c) Mask overlap (IoU)", f"All pairs (n = {len(pairs)})")

fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
print(f"Saved {OUT_PNG}")
