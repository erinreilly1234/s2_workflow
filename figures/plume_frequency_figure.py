"""
plume_frequency_figure.py
=========================

Makes ONE figure: two maps (S1 left, S2 right) showing, for every water pixel,
the percentage of scenes in which that pixel was flagged as plume.

    plume_areas_s1_s2.csv (step 07) ─► list of masks used ─┐
    S1 masks (resampled onto S2 grid) ─┐                    │
    S2 masks (step 06)                ─┴─► count per pixel ─┴─► % of scenes ─► maps
    water mask (step 05) ─► grey outside water
    outfall shapefile    ─► cyan triangles

The masks used are exactly the ones listed in plume_areas_s1_s2.csv, so this
figure always matches step 07 (same S2 rating filter). Masks with no plume
pixels still count as scenes (they add 0 to every pixel).

Output: fig_plume_frequency.png, saved next to this script (in figures/)

Needs: numpy, pandas, rasterio, geopandas, matplotlib
Run:   python plume_frequency_figure.py   (after 07_measure_plumes.py)
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling
import geopandas as gpd
import matplotlib
matplotlib.use("Agg")                     # write to file, no window
import matplotlib.pyplot as plt


# =============================================================================
# 1. FILE LOCATIONS
# =============================================================================
# TIDE server if /home/jovyan/s2 exists, otherwise the Mac external drive
# (set the TJ_ROOT environment variable if the drive is mounted elsewhere).
TIDE = Path("/home/jovyan/s2")
if TIDE.exists():
    AREAS_CSV   = TIDE / "07_figures" / "plume_areas_s1_s2.csv"
    S1_MASK_DIR = TIDE / "040_segoutput"
    S2_MASK_DIR = TIDE / "06_TURB_plumes"
    WATER_MASK  = TIDE / "water_mask_20240825.tif"
    OUTFALLS    = TIDE / "shapefiles" / "Outflow.shp"
else:
    TJ = Path(os.environ.get("TJ_ROOT", "/Volumes/External/TJ"))
    AREAS_CSV   = TJ / "010_optical" / "07_figures" / "plume_areas_s1_s2.csv"
    S1_MASK_DIR = TJ / "040_segoutput"
    S2_MASK_DIR = TJ / "010_optical" / "06_TURB_plumes"
    WATER_MASK  = TJ / "010_optical" / "water_mask_20240825.tif"
    OUTFALLS    = TJ / "015_shapefiles" / "outflow_PB_TJ" / "Outflow.shp"

OUT_PNG = Path(__file__).resolve().parent / "fig_plume_frequency.png"


# =============================================================================
# 2. SETTINGS
# =============================================================================
COLOUR_CAP_PERCENTILE = 99     # colour scale tops out at this percentile of non-zero values
COLORMAP = "magma_r"           # light = rarely plume, dark = often plume
LAND_COLOR = "#d9d9d9"         # grey for land / outside the water mask


# =============================================================================
# 3. HELPERS
# =============================================================================
def mask_on_grid(path, grid):
    """Plume mask (True/False) resampled onto the reference grid (nearest neighbour)."""
    out = np.zeros(grid["shape"], dtype="uint8")
    with rasterio.open(path) as src:
        reproject(rasterio.band(src, 1), out,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=grid["transform"], dst_crs=grid["crs"],
                  resampling=Resampling.nearest)
    return out == 1


def mask_paths(areas, sensor, folder):
    """Paths of this sensor's masks listed in step 07's CSV, found in `folder` by file name
    (so the CSV still works if it was written on another machine)."""
    names = areas.loc[areas["Sensor"] == sensor, "File"].map(lambda f: Path(f).name)
    return [folder / name for name in names]


def plume_frequency(paths, grid, water, resample):
    """% of scenes in which each water pixel is plume. Returns (map, number of scenes)."""
    count = np.zeros(grid["shape"], dtype=np.int32)
    n_scenes = 0
    for path in paths:
        if resample:                                   # S1: different grid -> resample
            plume = mask_on_grid(path, grid)
        else:                                          # S2: already on the grid
            with rasterio.open(path) as src:
                plume = src.read(1) == 1
            if plume.shape != grid["shape"]:
                print(f"  skipped (wrong size): {path.name}")
                continue
        count += plume
        n_scenes += 1
    percent = np.where(water, 100.0 * count / max(n_scenes, 1), np.nan)
    return percent, n_scenes


# =============================================================================
# 4. MAIN
# =============================================================================
def main():
    # --- reference grid and water mask (same grid as the S2 masks) ---
    with rasterio.open(WATER_MASK) as src:
        water = src.read(1) == 1
        grid = {"shape": water.shape, "transform": src.transform, "crs": src.crs}

    # --- frequency maps ---
    areas = pd.read_csv(AREAS_CSV)
    maps = {}
    for sensor, folder, resample in (("S1", S1_MASK_DIR, True), ("S2", S2_MASK_DIR, False)):
        maps[sensor] = plume_frequency(mask_paths(areas, sensor, folder), grid, water, resample)
        print(f"{sensor}: frequency from {maps[sensor][1]} scenes")

    # --- one colour scale for both maps, capped so a few hot pixels don't wash it out ---
    nonzero = np.concatenate([m[np.isfinite(m) & (m > 0)] for m, _ in maps.values()])
    vmax = np.nanpercentile(nonzero, COLOUR_CAP_PERCENTILE)

    # --- plot ---
    outfalls = gpd.read_file(OUTFALLS).to_crs(grid["crs"])
    cmap = plt.get_cmap(COLORMAP).copy()
    cmap.set_bad(LAND_COLOR)

    fig, axes = plt.subplots(1, 2, figsize=(13, 9))
    for ax, (sensor, (freq, n)) in zip(axes, maps.items()):
        image = ax.imshow(freq, cmap=cmap, vmin=0, vmax=vmax, interpolation="nearest")
        for point in outfalls.geometry:
            row, col = rasterio.transform.rowcol(grid["transform"], point.x, point.y)
            ax.plot(col, row, "^", color="cyan", mec="k", ms=8)
        ax.set_title(f"{sensor}: % of scenes flagged as plume (n = {n})")
        ax.axis("off")
    fig.colorbar(image, ax=axes, fraction=0.025, pad=0.02,
                 label=f"Plume frequency (% of scenes, colour scale capped at {vmax:.0f}%)")
    fig.suptitle("Plume frequency maps (grey = land / outside water mask)")

    fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
