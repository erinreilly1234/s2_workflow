"""
07_results_figures.py

Run after 06_turbidity_plumes.py. Makes the results figures from the S2 turbidity
plume masks and the S1 SAR plume masks:

  fig_area_flow_by_year.png   plume area through time (S1 and S2), one panel per
                              year, with the IBWC river flow (7-day rolling mean)
  fig_s1_s2_agreement.png     matched S1/S2 pairs: (a) area scatter with fit,
                              (b) plume centroid offset, (c) mask overlap (IoU)
  fig_plume_frequency.png     heat maps: % of scenes in which each pixel was
                              flagged as plume, S1 and S2 side by side (S1 is
                              resampled onto the S2 grid)
  plume_areas_s1_s2.csv       area of every mask (km²)
  matched_pairs_with_areas.csv  per-pair areas, centroid offset and IoU

S2 dates can be restricted to those rated cloud-free (Rating 3) in the ratings CSV.
Paths are for the jovyan server; copy the S1 masks and the three CSVs there first.
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from rasterio.warp import reproject, Resampling
from scipy.stats import linregress
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Paths ---
S2_DIR = Path("/home/jovyan/s2/06_TURB_plumes")            # *_TURB_mask.tif from step 06
S1_DIR = Path("/home/jovyan/s2/040_segoutput")             # S1 *_mask.tif
FLOW_CSV = Path("/home/jovyan/s2/ancillary/Total Flow.Preliminary Daily.US Million Gallons.csv")
RATINGS_CSV = Path("/home/jovyan/s2/all_s1_s2_dates_for_rating.csv")
PAIRS_CSV = Path("/home/jovyan/s2/ancillary/matched_pairs_reviewed.csv")
WATER_MASK = Path("/home/jovyan/s2/water_mask_20240825.tif")   # written by step 05
OUTFALLS = Path("/home/jovyan/s2/shapefiles/Outflow.shp")
OUT_DIR = Path("/home/jovyan/s2/07_figures")

# --- Settings ---
S2_RATING_3_ONLY = True      # keep only S2 dates rated cloud-free (Rating 3)
FLOW_ROLLING_DAYS = 7
DT = re.compile(r"(\d{8}T\d{6})")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def index_masks(folder, pattern):
    """{timestamp key 'YYYYMMDDTHHMMSS': path} for every mask file in folder."""
    out = {}
    for p in sorted(folder.glob(pattern)):
        m = DT.search(p.name)
        if m and not p.name.startswith("._"):
            out.setdefault(m.group(1), p)
    return out


def key_of(ts):
    return pd.to_datetime(ts, utc=True).strftime("%Y%m%dT%H%M%S")


def plume_area_km2(path):
    """Area (km²) of pixels == 1. Works for projected (m) and lat/lon grids."""
    with rasterio.open(path) as src:
        m = src.read(1) == 1
        t, crs = src.transform, src.crs
    if crs.is_geographic:                      # pixel area shrinks with latitude
        lat = t.f + t.e * (np.arange(m.shape[0]) + 0.5)
        row_area = abs(t.a) * 111320 * np.cos(np.radians(lat)) * abs(t.e) * 110540
        return float((m.sum(axis=1) * row_area).sum() / 1e6)
    return float(m.sum() * abs(t.a * t.e) / 1e6)


def mask_on_grid(path, ref):
    """Plume mask (bool) resampled onto the reference grid (nearest neighbour)."""
    out = np.zeros(ref["shape"], dtype="uint8")
    with rasterio.open(path) as src:
        reproject(rasterio.band(src, 1), out, src_transform=src.transform, src_crs=src.crs,
                  dst_transform=ref["transform"], dst_crs=ref["crs"],
                  resampling=Resampling.nearest)
    return out == 1


def load_flow(path):
    flow = pd.read_csv(path)
    ts = next(c for c in flow.columns if c.strip().startswith("Timestamp"))
    val = next(c for c in flow.columns if c.strip().startswith("Value"))
    flow = flow.rename(columns={ts: "Timestamp", val: "Flow"})
    flow["Timestamp"] = (pd.to_datetime(flow["Timestamp"], format="%m/%d/%y %H:%M")
                         .dt.tz_localize("Etc/GMT+8").dt.tz_convert("UTC"))
    flow = flow.sort_values("Timestamp")
    flow["Flow_roll"] = flow["Flow"].rolling(FLOW_ROLLING_DAYS, center=True, min_periods=1).mean()
    return flow


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    s2_index = index_masks(S2_DIR, "*_TURB_mask.tif")
    s1_index = index_masks(S1_DIR, "*_mask.tif")
    print(f"{len(s1_index)} S1 masks, {len(s2_index)} S2 masks")

    if S2_RATING_3_ONLY:
        r = pd.read_csv(RATINGS_CSV)
        good = set(r[(r["Sensor"] == "S2") & (r["Rating"] == 3)]["Timestamp"].map(key_of))
        s2_index = {k: v for k, v in s2_index.items() if k in good}
        print(f"  {len(s2_index)} S2 masks rated 3 (cloud-free)")

    # Reference grid = water mask (same grid as the S2 masks)
    with rasterio.open(WATER_MASK) as src:
        water = src.read(1) == 1
        ref = {"shape": water.shape, "transform": src.transform, "crs": src.crs}

    # ---- 1. Areas of every mask --------------------------------------------
    rows = []
    for sensor, idx in (("S1", s1_index), ("S2", s2_index)):
        for k, p in idx.items():
            rows.append({"Sensor": sensor, "Timestamp": pd.to_datetime(k, format="%Y%m%dT%H%M%S", utc=True),
                         "Area_km2": plume_area_km2(p), "File": str(p)})
    areas = pd.DataFrame(rows).sort_values("Timestamp")
    areas.to_csv(OUT_DIR / "plume_areas_s1_s2.csv", index=False)
    print("Saved plume_areas_s1_s2.csv")

    # ---- 2. Plume area vs river flow, by year ------------------------------
    flow = load_flow(FLOW_CSV)
    plotted = areas[areas["Area_km2"] > 0]               # dates with a plume
    years = sorted(plotted["Timestamp"].dt.year.unique())
    max_area = plotted["Area_km2"].max()
    max_flow = flow.loc[flow["Timestamp"].dt.year.isin(years), "Flow_roll"].max()

    fig, axes = plt.subplots(len(years), 1, figsize=(12, 4 * len(years)), squeeze=False)
    for ax, year in zip(axes[:, 0], years):
        for sensor, style in (("S1", dict(color="tab:blue", marker="o")),
                              ("S2", dict(color="tab:orange", marker="s"))):
            d = plotted[(plotted["Sensor"] == sensor) & (plotted["Timestamp"].dt.year == year)]
            label = "S1 plume area" if sensor == "S1" else \
                f"S2 plume area{' (rating 3 only)' if S2_RATING_3_ONLY else ''}"
            ax.plot(d["Timestamp"], d["Area_km2"], markersize=4, label=label, zorder=3, **style)
        ax.set_ylim(0, max_area * 1.1)
        ax.set_xlim(pd.Timestamp(f"{year}-01-01", tz="UTC"), pd.Timestamp(f"{year}-12-31", tz="UTC"))
        ax.set_ylabel("Plume area (km²)")
        ax.set_title(str(year), loc="left", fontweight="bold")
        ax2 = ax.twinx()
        f = flow[flow["Timestamp"].dt.year == year]
        ax2.plot(f["Timestamp"], f["Flow_roll"], color="gray", alpha=0.6, lw=1.2,
                 label=f"River flow ({FLOW_ROLLING_DAYS}-day rolling mean)", zorder=1)
        ax2.set_ylim(0, max_flow / 0.55)
        ax2.set_ylabel("River flow (Million Gal/day)")
    h1, l1 = axes[0, 0].get_legend_handles_labels()
    h2, l2 = fig.axes[1].get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.01), frameon=False)
    fig.suptitle("Plume area (S1 & S2) vs river flow, by year", y=1.03)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "fig_area_flow_by_year.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig_area_flow_by_year.png")

    # ---- 3. S1 vs S2 agreement on matched pairs ----------------------------
    pairs = pd.read_csv(PAIRS_CSV)
    out = []
    for _, row in pairs.iterrows():
        k1, k2 = key_of(row["S1_Time"]), key_of(row["S2_Time"])
        p1, p2 = s1_index.get(k1), s2_index.get(k2)
        rec = {"S1_Time": row["S1_Time"], "S2_Time": row["S2_Time"],
               "S2_Rating": row.get("S2_Rating"), "S1_File": p1, "S2_File": p2,
               "S1_Area_km2": np.nan, "S2_Area_km2": np.nan,
               "Centroid_Dist_m": np.nan, "IoU_percent": np.nan}
        if p1 and p2:
            m1 = mask_on_grid(p1, ref)
            with rasterio.open(p2) as s:
                m2 = s.read(1) == 1
            if m2.shape == m1.shape:
                px_km2 = abs(ref["transform"].a * ref["transform"].e) / 1e6
                rec["S1_Area_km2"], rec["S2_Area_km2"] = m1.sum() * px_km2, m2.sum() * px_km2
                union = (m1 | m2).sum()
                if union:
                    rec["IoU_percent"] = 100 * (m1 & m2).sum() / union
                if m1.any() and m2.any():
                    (r1, c1), (r2, c2) = np.argwhere(m1).mean(0), np.argwhere(m2).mean(0)
                    rec["Centroid_Dist_m"] = float(np.hypot(r1 - r2, c1 - c2) * abs(ref["transform"].a))
        out.append(rec)
    pairs_out = pd.DataFrame(out)
    pairs_out.to_csv(OUT_DIR / "matched_pairs_with_areas.csv", index=False)
    print("Saved matched_pairs_with_areas.csv")

    rated = pairs_out[pairs_out["S2_Rating"] == 3] if S2_RATING_3_ONLY else pairs_out
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(21, 6))
    fig.subplots_adjust(wspace=0.3)
    v = rated.dropna(subset=["S1_Area_km2", "S2_Area_km2"])
    a.scatter(v["S1_Area_km2"], v["S2_Area_km2"], color="steelblue", edgecolor="k", alpha=0.8)
    if len(v) >= 2:
        top = float(max(v["S1_Area_km2"].max(), v["S2_Area_km2"].max()))
        a.plot([0, top], [0, top], "--", color="gray", label="1:1 line")
        fit = linregress(v["S1_Area_km2"], v["S2_Area_km2"])
        xs = np.linspace(0, top, 200); ys = fit.intercept + fit.slope * xs
        keep = ys >= 0
        a.plot(xs[keep], ys[keep], color="firebrick", label="Best fit")
        a.set_xlim(0, top * 1.05); a.set_ylim(0, top * 1.05); a.set_aspect("equal")
        a.text(0.05, 0.95, f"$R^2$ = {fit.rvalue ** 2:.3f}\ny = {fit.slope:.3f}x + {fit.intercept:.2f}\nn = {len(v)}",
               transform=a.transAxes, va="top", fontsize=9,
               bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
        a.legend(fontsize=8, loc="lower right")
    a.set_xlabel("S1 plume area (km²)"); a.set_ylabel("S2 plume area (km²)")
    a.set_title("(a) S1 vs S2 plume area")
    # (b) and (c): distributions across all pairs (ECDF + individual pairs + summary stats)
    def ecdf_panel(ax, values, xlabel, title, unit, extra=""):
        v = np.sort(np.asarray(values, dtype=float))
        if len(v) == 0:
            ax.text(0.5, 0.5, "no pairs", ha="center", transform=ax.transAxes)
            return
        y = np.arange(1, len(v) + 1) / len(v)
        ax.step(np.r_[v[0], v], np.r_[0, y], where="post", color="0.3", lw=1.5, label="ECDF")
        ax.scatter(v, y, color="steelblue", edgecolor="k", s=40, zorder=3, label="individual pairs")
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.axvspan(q1, q3, color="steelblue", alpha=0.12, lw=0, label="interquartile range")
        ax.axvline(med, color="firebrick", ls="--", lw=1.2, label="median")
        ax.text(0.97, 0.05,
                f"n = {len(v)} pairs\nmedian = {med:.1f} {unit}\nIQR = {q1:.1f}–{q3:.1f} {unit}\n"
                f"mean = {v.mean():.1f} {unit}\nrange = {v.min():.1f}–{v.max():.1f} {unit}{extra}",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=9,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.9))
        ax.set_ylim(0, 1.02)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Cumulative fraction of pairs")
        ax.set_title(title)
        ax.legend(fontsize=8, loc="center right")

    d = rated["Centroid_Dist_m"].dropna() / 1000
    ecdf_panel(b, d, "Plume centroid distance, S1 vs S2 (km)", "(b) Centroid distance", "km")

    d = rated["IoU_percent"].dropna()
    ecdf_panel(c, d, "Plume mask overlap, S1 vs S2 (IoU, %)", "(c) Mask overlap (IoU)", "%",
               extra=f"\npairs with any overlap = {100 * (d > 0).mean():.0f}%" if len(d) else "")
    c.set_xlim(0, max(25, d.max() * 1.15) if len(d) else 100)

    fig.suptitle(f"S1 vs S2 plume agreement{' (S2 rating 3 only)' if S2_RATING_3_ONLY else ''}", y=1.02)
    fig.savefig(OUT_DIR / "fig_s1_s2_agreement.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig_s1_s2_agreement.png")

    # ---- 4. Plume frequency heat maps --------------------------------------
    freq = {}
    for sensor, idx in (("S1", s1_index), ("S2", s2_index)):
        count = np.zeros(ref["shape"], dtype=np.int32)
        n = 0
        for p in idx.values():
            if sensor == "S1":
                m = mask_on_grid(p, ref)
            else:
                with rasterio.open(p) as s:
                    m = s.read(1) == 1
                if m.shape != ref["shape"]:
                    continue
            count += m
            n += 1
        freq[sensor] = (np.where(water, 100.0 * count / max(n, 1), np.nan), n)
        print(f"  {sensor}: frequency from {n} scenes")

    vmax = np.nanpercentile(np.concatenate([f[np.isfinite(f) & (f > 0)] for f, _ in freq.values()]), 99)
    outfalls = gpd.read_file(OUTFALLS).to_crs(ref["crs"])
    fig, axes = plt.subplots(1, 2, figsize=(13, 9))
    cmap = plt.cm.magma_r.copy(); cmap.set_bad("#d9d9d9")
    for ax, (sensor, (f, n)) in zip(axes, freq.items()):
        im = ax.imshow(f, cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")
        for g in outfalls.geometry:
            r, cc = rasterio.transform.rowcol(ref["transform"], g.x, g.y)
            ax.plot(cc, r, "^", color="cyan", mec="k", ms=8)
        ax.set_title(f"{sensor}: % of scenes flagged as plume (n = {n})")
        ax.axis("off")
    fig.colorbar(im, ax=axes, fraction=0.025, pad=0.02,
                 label=f"Plume frequency (% of scenes, colour scale capped at {vmax:.0f}%)")
    fig.suptitle("Plume frequency maps (grey = land / outside water mask)")
    fig.savefig(OUT_DIR / "fig_plume_frequency.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved fig_plume_frequency.png")
    print("Done.")


if __name__ == "__main__":
    main()
