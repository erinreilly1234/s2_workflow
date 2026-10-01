"""
06_turbidity_plumes.py

Run after 05_compute_turbidity.py (uses the *_TURB.tif turbidity rasters, in FNU).

Same steps as 06_isegprob_ndti_926_2.py, applied to turbidity:
  1. smooth turbidity (low-pass filter that ignores NaN / land pixels)
  2. turbidity >= threshold marks turbid water
  3. morphological cleanup: opening, closing, fill small holes,
     remove small patches, dilation
  4. clip the cleaned mask to the step 05 water mask (smoothing, closing and
     dilation can otherwise push it onto land along the coast)
  5. distance from each pixel to the nearest outfall (meters)
  6. label clusters (4-connected); keep a cluster if any of its pixels lies
     within max_outfall_distance_m (1.5 km) of an outfall
Output per scene: plume mask (_TURB_mask.tif) and a 6-panel QC figure (qc_pngs/). Only S2 dates rated 3
(cloud-free) in RATINGS_CSV are processed when RATING_3_ONLY is True; other dates
are skipped entirely (no mask, no QC figure).

QC panel 1 colours pixels that have no turbidity value, using the red reflectance
from step 04 and the water mask from step 05:
  above the scale (magenta): turbidity > the colour-scale maximum, or red reflectance
                             >= C (too bright for the formula: cloud, sun glint, haze)
  below the scale (black):   red reflectance < 0 (very dark water, over-corrected)
  grey:                      land / outside the water mask / no data
"""

from pathlib import Path

import re

import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
from rasterio.features import rasterize
from scipy.ndimage import uniform_filter, distance_transform_edt, label, minimum
from rasterio.warp import reproject, Resampling
from skimage.morphology import (binary_opening, binary_closing, binary_dilation,
                                remove_small_objects, remove_small_holes, disk)
import matplotlib
matplotlib.use('Agg')  # write PNGs without needing a display
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- Parameters ---
input_dir = Path('/home/jovyan/s2/05_TURB')              # turbidity output of step 05
output_dir = Path('/home/jovyan/s2/06_TURB_plumes')
qc_dir = output_dir / 'qc_pngs'
outfall_shapefile = '/home/jovyan/s2/shapefiles/Outflow.shp'
reflectance_dir = Path('/home/jovyan/s2/03_mosaics_2022_2025_reflectance')   # step 04 output
water_mask_path = Path('/home/jovyan/s2/water_mask_20240825.tif')           # step 05 output
RATINGS_CSV = Path('/home/jovyan/s2/all_s1_s2_dates_for_rating.csv')
RATING_3_ONLY = True           # process only S2 dates rated 3 (cloud-free)
TURB_C = 0.2324                # same C as step 05 (red reflectance >= C has no turbidity)

threshold = 10.0               # turbidity (FNU) at or above this counts as turbid water
                               # (Nechad 2016 S2 coefficients from step 05; ~ the old 4 FNU)
max_outfall_distance_m = 1500.0  # keep a cluster if any pixel is within this distance of an outfall

# Spatial settings in meters / square meters (converted to pixels below)
smooth_window_m = 50.0         # low-pass filter window width
open_radius_m = 30.0           # opening: removes specks/thin features
close_radius_m = 50.0          # closing: bridges small gaps
min_hole_area_m2 = 80000.0     # fill holes smaller than this
min_object_area_m2 = 30000.0   # drop patches smaller than this
dilate_radius_m = 30.0         # final dilation buffer


with rasterio.open(water_mask_path) as _src:
    WATER = _src.read(1) == 1
    WATER_TRANSFORM, WATER_CRS = _src.transform, _src.crs


def water_on_grid(shape, transform, crs):
    """Step 05 water mask on a scene's pixel grid (same logic as step 05: resampled,
    nearest neighbour, when a scene's grid differs; outside the reference = not water)."""
    if shape == WATER.shape and transform == WATER_TRANSFORM:
        return WATER
    out = np.zeros(shape, dtype='uint8')
    reproject(WATER.astype('uint8'), out, src_transform=WATER_TRANSFORM, src_crs=WATER_CRS,
              dst_transform=transform, dst_crs=crs, resampling=Resampling.nearest)
    return out.astype(bool)


def process_image(image_path, mask_path, qc_path):
    # 1. Read the turbidity raster (land/no data = NaN)
    with rasterio.open(image_path) as src:
        turb = src.read(1).astype(np.float32)
        profile = src.profile
        transform = src.transform
        crs = src.crs

    px = abs(transform.a)                      # pixel size in meters (UTM)
    m_to_px = lambda meters: max(1, round(meters / px))
    m2_to_px = lambda m2: max(1, round(m2 / (px * px)))

    # 2. Low-pass filter, ignoring NaN pixels
    win = m_to_px(smooth_window_m)
    if win % 2 == 0:
        win += 1                               # odd window stays centered
    valid = ~np.isnan(turb)
    num = uniform_filter(np.where(valid, turb, 0.0), size=win, mode='nearest')
    den = uniform_filter(valid.astype(np.float32), size=win, mode='nearest')
    smooth = np.where(den > 0, num / np.maximum(den, 1e-6), np.nan)

    # 3. Threshold & morphological cleanup
    signal = smooth >= threshold               # NaN compares False
    clean = binary_opening(signal, disk(m_to_px(open_radius_m)))
    clean = binary_closing(clean, disk(m_to_px(close_radius_m)))
    clean = remove_small_holes(clean, area_threshold=m2_to_px(min_hole_area_m2))
    clean = remove_small_objects(clean, min_size=m2_to_px(min_object_area_m2))
    mask = binary_dilation(clean, disk(m_to_px(dilate_radius_m)))

    # 3b. Clip to the water mask from step 05 (on this scene's grid)
    water = water_on_grid(turb.shape, transform, crs)
    mask &= water

    # 4. Distance (meters) from every pixel to the nearest outfall
    outfalls = gpd.read_file(outfall_shapefile).to_crs(crs)
    outfall_r = rasterize([(g, 1) for g in outfalls.geometry], out_shape=turb.shape,
                          transform=transform, fill=0, dtype='uint8')
    dist_m = distance_transform_edt(outfall_r == 0) * px

    if not outfall_r.any():
        raise ValueError(f'No outfall from {outfall_shapefile} falls inside {image_path.name}')

    # 5. Label clusters; keep those with any pixel within max_outfall_distance_m of an outfall
    labels, n_clusters = label(mask)
    mask_selected = np.zeros(turb.shape, dtype=bool)
    n_accepted = 0
    if n_clusters > 0:
        ids = np.arange(1, n_clusters + 1)
        cluster_dist = np.asarray(minimum(dist_m, labels, index=ids), dtype=np.float32)
        accepted = ids[cluster_dist <= max_outfall_distance_m]
        n_accepted = len(accepted)
        mask_selected = np.isin(labels, accepted)

    # 6. Write the plume mask
    out_profile = profile.copy()
    out_profile.update(count=1, dtype='uint8', nodata=0, compress='deflate')
    with rasterio.open(mask_path, 'w', **out_profile) as dst:
        dst.write(mask_selected.astype('uint8'), 1)

    area_km2 = mask_selected.sum() * px * px / 1e6
    print(f'  plume area: {area_km2:.2f} km2 ({n_accepted} of {n_clusters} clusters kept)')

    # 8. QC figure (6 panels, same layout as the NDTI scripts)
    if qc_path is None:
        return
    ys, xs = np.where(outfall_r == 1)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # Panel 1: turbidity, with pixels that have no value shown as above / below the scale
    vmax = 3 * threshold
    disp = turb.copy()
    stem = image_path.name.replace('_TURB.tif', '')
    refl_path = reflectance_dir / f'{stem}_reflectance.tif'
    below = above_sat = np.zeros(turb.shape, dtype=bool)
    if refl_path.exists() and water.any():
        with rasterio.open(refl_path) as src:
            red = src.read(4)
        if red.shape == turb.shape:
            below = water & (red < 0)                 # very dark / over-corrected water
            above_sat = water & (red >= TURB_C)       # too bright: cloud, glint, haze
    disp[below] = -1.0                                # shown in the 'under' colour
    disp[above_sat] = vmax * 10                       # shown in the 'over' colour
    n_water = max(water.sum(), 1)
    pct_over = 100 * (water & (np.nan_to_num(turb) > vmax)).sum() / n_water
    pct_sat = 100 * above_sat.sum() / n_water
    pct_below = 100 * below.sum() / n_water

    cmap = plt.cm.viridis.copy()
    cmap.set_under('black')
    cmap.set_over('#ff3bd3')
    cmap.set_bad('#d9d9d9')
    im1 = axes[0, 0].imshow(disp, cmap=cmap, vmin=0, vmax=vmax, interpolation='nearest')
    fig.colorbar(im1, ax=axes[0, 0], fraction=0.046, extend='both', label='Turbidity (FNU)')
    axes[0, 0].legend(handles=[
        Patch(color='#ff3bd3', label=f'above scale: > {vmax:g} FNU ({pct_over:.1f}% of water)\n'
                                     f'or too bright, cloud/glint/haze ({pct_sat:.1f}%)'),
        Patch(color='black', label=f'below scale: red reflectance < 0 ({pct_below:.1f}%)'),
        Patch(color='#d9d9d9', label='land / no data'),
    ], loc='lower left', fontsize=7, framealpha=0.9)
    axes[0, 0].set_title(f'1. Turbidity (FNU, colour scale 0-{vmax:g})')
    axes[0, 1].imshow(signal, cmap='gray')
    axes[0, 1].set_title(f'2. Turbidity >= {threshold:g} FNU\n{signal.sum():,} px')
    axes[0, 2].imshow(mask, cmap='gray')
    axes[0, 2].set_title(f'3. After morphology\n{mask.sum():,} px')
    im = axes[1, 0].imshow(dist_m / 1000, cmap='inferno_r', vmin=0, vmax=3 * max_outfall_distance_m / 1000)
    axes[1, 0].contour(dist_m, levels=[max_outfall_distance_m], colors='cyan', linewidths=1)
    axes[1, 0].scatter(xs, ys, c='cyan', s=15, edgecolor='k')
    axes[1, 0].set_title(f'4. Distance to outfall (km); cyan = {max_outfall_distance_m / 1000:g} km')
    fig.colorbar(im, ax=axes[1, 0], fraction=0.046)
    axes[1, 1].imshow(labels, cmap='nipy_spectral', interpolation='nearest')
    axes[1, 1].set_title(f'5. Clusters (n={n_clusters}, {n_accepted} accepted)')
    axes[1, 2].imshow(turb, cmap='gray', vmin=0, vmax=3 * threshold)
    axes[1, 2].imshow(np.ma.masked_where(~mask_selected, mask_selected),
                      cmap='autumn', alpha=0.6, interpolation='nearest')
    axes[1, 2].scatter(xs, ys, c='cyan', s=15, edgecolor='k')
    axes[1, 2].set_title(f'6. Final plume mask\n{area_km2:.2f} km2')
    for ax in axes.ravel():
        ax.axis('off')
    fig.suptitle(image_path.name, fontsize=13)
    fig.tight_layout()
    fig.savefig(qc_path, dpi=130, bbox_inches='tight')
    plt.close(fig)


def rated_3_keys():
    """Timestamps (YYYYMMDDTHHMMSS, UTC) of S2 dates rated 3 in RATINGS_CSV."""
    r = pd.read_csv(RATINGS_CSV)
    r = r[(r['Sensor'] == 'S2') & (r['Rating'] == 3)]
    return set(pd.to_datetime(r['Timestamp'], utc=True).dt.strftime('%Y%m%dT%H%M%S'))


def main():
    qc_dir.mkdir(parents=True, exist_ok=True)
    good = rated_3_keys() if RATING_3_ONLY else None
    if good is not None:
        print(f'Processing only the {len(good)} S2 dates rated 3')
    for src_path in sorted(input_dir.glob('*_TURB.tif')):
        if src_path.name.startswith('._'):     # macOS metadata files on external drives
            continue
        m = re.search(r'(\d{8}T\d{6})', src_path.name)
        if good is not None and not (m and m.group(1) in good):
            continue                           # not rated 3: skip this date
        mask_path = output_dir / src_path.name.replace('_TURB.tif', '_TURB_mask.tif')
        qc_path = qc_dir / src_path.name.replace('_TURB.tif', '_TURB_qc.png')

        if mask_path.exists() and qc_path.exists():
            print(f'Skipping {src_path.name} (already done)')
            continue

        print(f'Processing {src_path.name}')
        process_image(src_path, mask_path, qc_path)

    print('Done.')


if __name__ == '__main__':
    main()
