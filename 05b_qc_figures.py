"""
QC figures for the Sentinel-2 reflectance and NDTI mosaics.

Run after 05_compute_ndti_ndci.py (independent of the 06 segmentation step). For each date writes 06_QC/QC_<date>.png with:
  - true colour (B04/B03/B02), fixed linear stretch 0-0.3 reflectance, gamma 1/1.6
    (same stretch on every date, so haze/cloud shows up as brighter)
  - spectral profiles (B01-B8A) at POINTS, with 3x3 window min-max shading
  - NDTI map, colour stretch = that date's 2nd-98th percentile
  - NDTI histogram with p2 / median / p98 and the point values
Also writes QC_overview_all_dates.png and QC_stats.csv (NDTI stats + point spectra).

Band order in the mosaics: B01 B02 B03 B04 B05 B06 B07 B08 B8A B11 B12.
"""
from pathlib import Path
import numpy as np, rasterio, csv, gc
from rasterio.warp import transform as warp_transform
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
HALO = [pe.withStroke(linewidth=2.5, foreground='white')]

refl_dir = Path("/home/jovyan/s2/03_mosaics_2022_2025_reflectance")
ndti_dir = Path("/home/jovyan/s2/05_NDTI")
out = Path("/home/jovyan/s2/06_QC")
out.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#888", "axes.labelcolor": "#333", "xtick.color": "#555", "ytick.color": "#555"})
POINTS = [("P1", -117.141724, 32.556363, "#e03131", "o"), ("P2", -117.220688, 32.552601, "#f08c00", "s")]
BANDS = ["B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B8A", "B11", "B12"]
WL = np.array([443, 492, 560, 665, 704, 740, 783, 833, 865, 1614, 2202])
RGB_MAX, GAMMA = 0.30, 1/1.6
NDTI_CMAP = "viridis"
WATER_MASK_LABEL = "2024-08-25"   # keep in sync with WATER_MASK_DATE in step 05
rows, thumbs = [], []

def date_of(name): return name.split("_")[1][:8]

for rp in sorted(p for p in refl_dir.glob("*.tif") if not p.name.startswith("._")):
    d = date_of(rp.name)
    with rasterio.open(rp) as s:
        xs, ys = warp_transform("EPSG:4326", s.crs, [p[1] for p in POINTS], [p[2] for p in POINTS])
        pts = []
        for (name, lon, lat, col, mk), x, y in zip(POINTS, xs, ys):
            r, c = s.index(x, y)
            if 0 <= r < s.height and 0 <= c < s.width:
                win = s.read(window=((r-1, r+2), (c-1, c+2)), boundless=True, fill_value=np.nan)
                pts.append((name, lon, lat, col, mk, r, c, win))
        rgb = np.dstack([s.read(4), s.read(3), s.read(2)])
    with rasterio.open(ndti_dir / rp.name.replace("_reflectance.tif", "_NDTI.tif")) as s:
        ndti = s.read(1)
    valid = np.all(np.isfinite(rgb), axis=2)
    img = (np.clip(np.nan_to_num(rgb) / RGB_MAX, 0, 1) ** GAMMA).astype(np.float32); del rgb
    img = np.dstack([img, valid.astype(np.float32)])
    nv = ndti[np.isfinite(ndti)]
    q = np.percentile(nv, [2, 50, 98]) if nv.size else np.array([np.nan]*3)

    stats = {"date": d, "file": rp.name, "ndti_n": int(nv.size),
             "ndti_p2": round(float(q[0]), 4), "ndti_median": round(float(q[1]), 4), "ndti_p98": round(float(q[2]), 4),
             "ndti_mean": round(float(nv.mean()), 4) if nv.size else "", "ndti_sd": round(float(nv.std()), 4) if nv.size else "",
             }
    for name, lon, lat, col, mk, r, c, win in pts:
        stats[f"{name}_row"], stats[f"{name}_col"] = r, c
        stats[f"{name}_ndti"] = round(float(ndti[r, c]), 4)
        for b, v in zip(BANDS, win[:, 1, 1]): stats[f"{name}_{b}"] = round(float(v), 5)

    def mark(a, scale=1, ms=9, label=True):
        for name, lon, lat, col, mk, r, c, win in pts:
            a.plot(c/scale, r/scale, marker=mk, ms=ms, mfc="none", mec=col, mew=2 if ms > 5 else 1, path_effects=[pe.withStroke(linewidth=4, foreground="white")])
            if label: a.annotate(name, (c/scale, r/scale), xytext=(-8, 0), textcoords="offset points", ha="right", va="center",
                                 color=col, fontsize=8, fontweight="bold", path_effects=HALO)

    fig, ax = plt.subplots(2, 2, figsize=(12, 11), gridspec_kw={"width_ratios": [1, 1.15]})
    fig.suptitle(f"Sentinel-2 QC  ·  {d[:4]}-{d[4:6]}-{d[6:]}  ·  {rp.name}", fontsize=11, x=0.02, ha="left")
    a = ax[0, 0]; a.imshow(img); a.set_title(f"True colour (B04/B03/B02), linear 0–{RGB_MAX} refl., γ={GAMMA:.2f}", loc="left"); a.axis("off")
    mark(a)

    a = ax[0, 1]
    NP = 9  # plot B01-B8A only (skip B11, B12)
    for name, lon, lat, col, mk, r, c, win in pts:
        flat = win[:NP].reshape(NP, -1)
        a.fill_between(WL[:NP], np.nanmin(flat, 1), np.nanmax(flat, 1), color=col, alpha=0.13, lw=0)
        a.plot(WL[:NP], win[:NP, 1, 1], "-", marker=mk, color=col, lw=2, ms=6, label=f"{name}  {lat:.5f}, {lon:.5f}")
    a.set_xticks(WL[:NP]); a.set_xticklabels([f"{b}\n{w}" for b, w in zip(BANDS[:NP], WL[:NP])], fontsize=7.5)
    a.axhline(0, color="#999", lw=0.8, ls="--")
    a.legend(frameon=False, fontsize=8, loc="upper right", title="pixel value; shading = 3×3 min–max", title_fontsize=7.5)
    if not pts: a.text(0.5, 0.5, "Points outside mosaic", ha="center", transform=a.transAxes)
    a.set_xlabel("Band / wavelength (nm)"); a.set_ylabel("BOA reflectance")
    a.set_title("Spectral profiles", loc="left")

    a = ax[1, 0]
    im = a.imshow(ndti, cmap=NDTI_CMAP, vmin=q[0], vmax=q[2], interpolation="nearest")
    mark(a)
    a.set_title(f"NDTI, fixed water mask ({WATER_MASK_LABEL}) · stretch p2–p98", loc="left"); a.axis("off")
    cb = fig.colorbar(im, ax=a, fraction=0.035, pad=0.02); cb.set_label(f"NDTI  [{q[0]:.3f} to {q[2]:.3f}]")

    a = ax[1, 1]
    if nv.size:
        lo, hi = np.percentile(nv, [0.5, 99.5]); pad = 0.15*(hi-lo)
        a.hist(nv, bins=np.linspace(lo-pad, hi+pad, 241), color="#5b8a8a", edgecolor="none")
        a.axvspan(q[0], q[2], color="#5b8a8a", alpha=0.08, lw=0)
        for x, l in zip(q, ["p2", "median", "p98"]):
            a.axvline(x, color="#333" if l == "median" else "#999", lw=1, ls="-" if l == "median" else ":")
        for name, lon, lat, col, mk, r, c, win in pts:
            if np.isfinite(ndti[r, c]): a.axvline(ndti[r, c], color=col, lw=1.5, label=f"{name} {ndti[r, c]:.3f}")
        if pts: a.legend(frameon=False, loc="upper left", fontsize=8)
        a.text(0.98, 0.95, f"n = {nv.size:,}\nmedian {q[1]:.3f}\np2 {q[0]:.3f}  p98 {q[2]:.3f}\nmean {nv.mean():.3f}  sd {nv.std():.3f}",
               transform=a.transAxes, ha="right", va="top", fontsize=8.5, color="#333")
    a.set_xlabel("NDTI (x-range p0.5–p99.5)"); a.set_ylabel("pixel count"); a.set_title("NDTI distribution (water pixels)", loc="left")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / f"QC_{d}.png", dpi=110); plt.close(fig)
    rows.append(stats); thumbs.append((d, img[::3, ::3].copy(), ndti[::3, ::3].copy(), q[0], q[2], [(p[3], p[4], p[5]//3, p[6]//3) for p in pts]))
    print(d, stats.get("P1_ndti"), stats.get("P2_ndti"), flush=True)
    del img, ndti, pts; gc.collect()

keys = list(dict.fromkeys(k for r_ in rows for k in r_))
with open(out/"QC_stats.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)

n = len(thumbs); fig, ax = plt.subplots(2, n, figsize=(1.9*n, 6.4), gridspec_kw={"wspace": 0.05})
for j, (d, img, nd, lo, hi, pp) in enumerate(thumbs):
    ax[0, j].imshow(img); ax[0, j].set_title(f"{d[:4]}-{d[4:6]}-{d[6:]}", fontsize=8)
    ax[1, j].imshow(nd, cmap=NDTI_CMAP, vmin=lo, vmax=hi)
    ax[1, j].set_title(f"{lo:.2f} to {hi:.2f}", fontsize=7, color="#555")
    for a in ax[:, j]:
        for col, mk, rr, cc in pp: a.plot(cc, rr, marker=mk, ms=4, mfc="none", mec=col, mew=1)
        a.axis("off")
fig.suptitle("All dates — true colour (top, same stretch) and NDTI (bottom, each date stretched p2–p98; range shown above each)", x=0.01, ha="left", fontsize=10)
fig.savefig(out/"QC_overview_all_dates.png", dpi=120, bbox_inches="tight")
print("overview done")
