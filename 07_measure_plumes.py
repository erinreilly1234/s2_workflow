"""
07_measure_plumes.py
====================

Step 07 of the workflow. Measures every plume mask (S1 and S2) and compares the
S1 and S2 masks of matched same-day pairs. This step only calculates numbers;
the figures are made by the scripts in 010_optical/08_figures/ (see the README).

    S2 plume masks (step 06) ─┐
    S1 plume masks (SAR)     ─┼─► plume_areas_s1_s2.csv        one row per mask
    S2 rating spreadsheet    ─┤
    matched S1/S2 pairs      ─┴─► matched_pairs_with_areas.csv one row per pair

Outputs (written to 010_optical/07_plume_measurements/)
------------------------------------------------------
plume_areas_s1_s2.csv
    Sensor, Timestamp (UTC), Area_km2, File
    Area of every mask. Area_km2 = 0 means the mask has no plume pixels.

matched_pairs_with_areas.csv
    S1_Time, S2_Time, S2_Rating, S1_File, S2_File,
    S1_Area_km2, S2_Area_km2, Centroid_Dist_m, IoU_percent
    For each pair in matched_pairs_reviewed.csv: both masks on the same
    (S2) grid, their areas, the distance between plume centres and the
    mask overlap (intersection over union). Blank if a mask is missing.

Which S2 masks are used
-----------------------
Only S2 dates rated 3 (cloud-free) in all_s1_s2_dates_for_rating.csv, when
S2_RATING_3_ONLY = True. All S1 masks are used (radar sees through cloud).

Needs: numpy, pandas, rasterio
Run:   python 07_measure_plumes.py
"""

import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling


# =============================================================================
# 1. FILE LOCATIONS
# =============================================================================
# TIDE server if /home/jovyan/s2 exists, otherwise the Mac external drive
# (set the TJ_ROOT environment variable if the drive is mounted elsewhere).
TIDE = Path("/home/jovyan/s2")
if TIDE.exists():
    S2_MASK_DIR = TIDE / "06_TURB_plumes"
    S1_MASK_DIR = TIDE / "040_segoutput"
    RATINGS_CSV = TIDE / "all_s1_s2_dates_for_rating.csv"
    PAIRS_CSV   = TIDE / "ancillary" / "matched_pairs_reviewed.csv"
    WATER_MASK  = TIDE / "water_mask_20240825.tif"
    OUT_DIR     = TIDE / "07_plume_measurements"
else:
    TJ = Path(os.environ.get("TJ_ROOT", "/Volumes/External/TJ"))
    S2_MASK_DIR = TJ / "010_optical" / "06_TURB_plumes"        # *_TURB_mask.tif from step 06
    S1_MASK_DIR = TJ / "040_segoutput"                         # S1 *_mask.tif
    RATINGS_CSV = TJ / "__spreadsheets" / "all_s1_s2_dates_for_rating.csv"
    PAIRS_CSV   = TJ / "__spreadsheets" / "matched_pairs_reviewed.csv"
    WATER_MASK  = TJ / "010_optical" / "water_mask_20240825.tif"   # step 05; defines the S2 grid
    OUT_DIR     = TJ / "010_optical" / "07_plume_measurements"

AREAS_OUT = OUT_DIR / "plume_areas_s1_s2.csv"
PAIRS_OUT = OUT_DIR / "matched_pairs_with_areas.csv"


# =============================================================================
# 2. SETTINGS
# =============================================================================
S2_RATING_3_ONLY = True                  # keep only S2 dates rated cloud-free (Rating 3)
TIME_IN_NAME = re.compile(r"(\d{8}T\d{6})")   # e.g. 20220101T015007 in a file name


# =============================================================================
# 3. FIND THE MASK FILES
# =============================================================================
def find_masks(folder, pattern):
    """{'YYYYMMDDTHHMMSS': path} for every mask in folder (first file per time)."""
    masks = {}
    for path in sorted(folder.glob(pattern)):
        if path.name.startswith("._"):          # macOS metadata files on external drives
            continue
        match = TIME_IN_NAME.search(path.name)
        if match:
            masks.setdefault(match.group(1), path)
    return masks


def time_key(timestamp):
    """Any timestamp -> 'YYYYMMDDTHHMMSS' (UTC), the key used by find_masks."""
    return pd.to_datetime(timestamp, utc=True).strftime("%Y%m%dT%H%M%S")


def cloud_free_s2_keys():
    """Time keys of S2 scenes rated 3 (cloud-free) in the ratings spreadsheet."""
    ratings = pd.read_csv(RATINGS_CSV)
    rated_3 = ratings[(ratings["Sensor"] == "S2") & (ratings["Rating"] == 3)]
    return set(rated_3["Timestamp"].map(time_key))


# =============================================================================
# 4. MEASURE A MASK
# =============================================================================
def plume_area_km2(path):
    """Area (km²) of the plume pixels (value 1). Works for metre and lat/lon grids."""
    with rasterio.open(path) as src:
        plume = src.read(1) == 1
        t, crs = src.transform, src.crs
    if crs.is_geographic:
        # lat/lon grid: pixel width in metres shrinks with latitude, so do it row by row
        lat = t.f + t.e * (np.arange(plume.shape[0]) + 0.5)
        row_pixel_area_m2 = abs(t.a) * 111320 * np.cos(np.radians(lat)) * abs(t.e) * 110540
        return float((plume.sum(axis=1) * row_pixel_area_m2).sum() / 1e6)
    return float(plume.sum() * abs(t.a * t.e) / 1e6)


def mask_on_grid(path, grid):
    """Plume mask (True/False) resampled onto the reference grid (nearest neighbour)."""
    out = np.zeros(grid["shape"], dtype="uint8")
    with rasterio.open(path) as src:
        reproject(rasterio.band(src, 1), out,
                  src_transform=src.transform, src_crs=src.crs,
                  dst_transform=grid["transform"], dst_crs=grid["crs"],
                  resampling=Resampling.nearest)
    return out == 1


# =============================================================================
# 5. COMPARE THE TWO MASKS OF A PAIR
# =============================================================================
def compare_masks(s1_mask, s2_mask, pixel_size_m):
    """Areas, centre-to-centre distance and overlap (IoU) of two masks on one grid."""
    pixel_area_km2 = pixel_size_m ** 2 / 1e6
    result = {"S1_Area_km2": s1_mask.sum() * pixel_area_km2,
              "S2_Area_km2": s2_mask.sum() * pixel_area_km2,
              "Centroid_Dist_m": np.nan,
              "IoU_percent": np.nan}

    union = (s1_mask | s2_mask).sum()
    if union:                                   # at least one sensor found a plume
        result["IoU_percent"] = 100 * (s1_mask & s2_mask).sum() / union

    if s1_mask.any() and s2_mask.any():         # both found a plume
        r1, c1 = np.argwhere(s1_mask).mean(axis=0)
        r2, c2 = np.argwhere(s2_mask).mean(axis=0)
        result["Centroid_Dist_m"] = float(np.hypot(r1 - r2, c1 - c2) * pixel_size_m)
    return result


# =============================================================================
# 6. MAIN
# =============================================================================
def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # --- which masks exist ---
    s1_masks = find_masks(S1_MASK_DIR, "*_mask.tif")
    s2_masks = find_masks(S2_MASK_DIR, "*_TURB_mask.tif")
    print(f"Found {len(s1_masks)} S1 masks and {len(s2_masks)} S2 masks")
    if S2_RATING_3_ONLY:
        keep = cloud_free_s2_keys()
        s2_masks = {k: p for k, p in s2_masks.items() if k in keep}
        print(f"  kept {len(s2_masks)} S2 masks rated 3 (cloud-free)")

    # --- A. area of every mask -> plume_areas_s1_s2.csv ---
    rows = []
    for sensor, masks in (("S1", s1_masks), ("S2", s2_masks)):
        for key, path in masks.items():
            rows.append({"Sensor": sensor,
                         "Timestamp": pd.to_datetime(key, format="%Y%m%dT%H%M%S", utc=True),
                         "Area_km2": plume_area_km2(path),
                         "File": str(path)})
    areas = pd.DataFrame(rows).sort_values("Timestamp")
    areas.to_csv(AREAS_OUT, index=False)
    print(f"Saved {AREAS_OUT.name} ({len(areas)} masks)")

    # --- B. matched S1/S2 pairs -> matched_pairs_with_areas.csv ---
    # Reference grid = the water mask (same grid as the S2 masks); S1 is resampled onto it.
    with rasterio.open(WATER_MASK) as src:
        grid = {"shape": (src.height, src.width), "transform": src.transform, "crs": src.crs}
    pixel_size_m = abs(grid["transform"].a)

    rows = []
    for _, pair in pd.read_csv(PAIRS_CSV).iterrows():
        s1_path = s1_masks.get(time_key(pair["S1_Time"]))
        s2_path = s2_masks.get(time_key(pair["S2_Time"]))   # None if not rated 3 or no mask
        row = {"S1_Time": pair["S1_Time"], "S2_Time": pair["S2_Time"],
               "S2_Rating": pair.get("S2_Rating"),
               "S1_File": s1_path, "S2_File": s2_path,
               "S1_Area_km2": np.nan, "S2_Area_km2": np.nan,
               "Centroid_Dist_m": np.nan, "IoU_percent": np.nan}
        if s1_path and s2_path:
            s1_mask = mask_on_grid(s1_path, grid)
            with rasterio.open(s2_path) as src:
                s2_mask = src.read(1) == 1
            if s2_mask.shape == s1_mask.shape:
                row.update(compare_masks(s1_mask, s2_mask, pixel_size_m))
        rows.append(row)
    pairs = pd.DataFrame(rows)
    pairs.to_csv(PAIRS_OUT, index=False)
    n_both = pairs["S1_Area_km2"].notna().sum()
    print(f"Saved {PAIRS_OUT.name} ({len(pairs)} pairs, {n_both} with both masks)")


if __name__ == "__main__":
    main()
