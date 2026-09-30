"""
06_turbidity_plumes.py

Run after 05_compute_ndti_ndci.py (uses the *_TURB.tif turbidity rasters, in FNU).

Same steps as 06_isegprob_ndti_926_2.py, applied to turbidity:
  1. smooth turbidity (low-pass filter that ignores NaN / land pixels)
  2. turbidity >= threshold marks turbid water
  3. morphological cleanup: opening, closing, fill small holes,
     remove small patches, dilation
  4. distance from each pixel to the nearest outfall (meters)
  5. outfall probability = exp(-distance / decay_scale)
  6. label clusters; each cluster gets the highest probability inside it;
     keep clusters with probability >= cluster_prob_threshold
Outputs per scene: plume mask (_TURB_mask.tif), cluster probability map
(_TURB_prob.tif) and a 6-panel QC figure (qc_pngs/).
"""

from pathlib import Path

import numpy as np
import rasterio
import geopandas as gpd
from rasterio.features import rasterize
from scipy.ndimage import uniform_filter, distance_transform_edt, label, maximum
from skimage.morphology import (binary_opening, binary_closing, binary_dilation,
                                remove_small_objects, remove_small_holes, disk)
import matplotlib
matplotlib.use('Agg')  # write PNGs without needing a display
import matplotlib.pyplot as plt

# --- Parameters ---
input_dir = Path('/home/jovyan/s2/05_TURB')              # turbidity output of step 05
output_dir = Path('/home/jovyan/s2/06_TURB_plumes')
qc_dir = output_dir / 'qc_pngs'
outfall_shapefile = '/home/jovyan/s2/shapefiles/Outflow.shp'

threshold = 10.0               # turbidity (FNU) at or above this counts as turbid water
                               # (Nechad 2016 S2 coefficients from step 05; ~ the old 4 FNU)
decay_scale = 4000.0           # distance decay scale (meters)
cluster_prob_threshold = 0.7   # minimum probability for cluster acceptance

# Spatial settings in meters / square meters (converted to pixels below)
smooth_window_m = 50.0         # low-pass filter window width
open_radius_m = 30.0           # opening: removes specks/thin features
close_radius_m = 50.0          # closing: bridges small gaps
min_hole_area_m2 = 80000.0     # fill holes smaller than this
min_object_area_m2 = 30000.0   # drop patches smaller than this
dilate_radius_m = 30.0         # final dilation buffer


def process_image(image_path, mask_path, prob_path, qc_path):
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

    # 4. Distance (meters) from every pixel to the nearest outfall
    outfalls = gpd.read_file(outfall_shapefile).to_crs(crs)
    outfall_r = rasterize([(g, 1) for g in outfalls.geometry], out_shape=turb.shape,
                          transform=transform, fill=0, dtype='uint8')
    dist_m = distance_transform_edt(outfall_r == 0) * px

    # 5. Exponential decay probability
    prob_full = np.exp(-dist_m / decay_scale).astype(np.float32)

    # 6. Label clusters, give each its highest probability, keep the likely ones
    labels, n_clusters = label(mask)
    prob_map = np.zeros(turb.shape, dtype=np.float32)
    mask_selected = np.zeros(turb.shape, dtype=bool)
    n_accepted = 0
    if n_clusters > 0:
        ids = np.arange(1, n_clusters + 1)
        cluster_prob = np.asarray(maximum(prob_full, labels, index=ids), dtype=np.float32)
        prob_map = np.concatenate([[0], cluster_prob])[labels].astype(np.float32)
        accepted = ids[cluster_prob >= cluster_prob_threshold]
        n_accepted = len(accepted)
        mask_selected = np.isin(labels, accepted)

    # 7. Write the plume mask and the cluster probability map
    out_profile = profile.copy()
    out_profile.update(count=1, dtype='uint8', nodata=0, compress='deflate')
    with rasterio.open(mask_path, 'w', **out_profile) as dst:
        dst.write(mask_selected.astype('uint8'), 1)
    out_profile.update(dtype='float32')
    with rasterio.open(prob_path, 'w', **out_profile) as dst:
        dst.write(prob_map, 1)

    area_km2 = mask_selected.sum() * px * px / 1e6
    print(f'  plume area: {area_km2:.2f} km2 ({n_accepted} of {n_clusters} clusters kept)')

    # 8. QC figure (6 panels, same layout as the NDTI scripts)
    ys, xs = np.where(outfall_r == 1)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    axes[0, 0].imshow(turb, cmap='viridis', vmin=0, vmax=3 * threshold)
    axes[0, 0].set_title(f'1. Turbidity (FNU, colour scale 0-{3 * threshold:g})')
    axes[0, 1].imshow(signal, cmap='gray')
    axes[0, 1].set_title(f'2. Turbidity >= {threshold:g} FNU\n{signal.sum():,} px')
    axes[0, 2].imshow(mask, cmap='gray')
    axes[0, 2].set_title(f'3. After morphology\n{mask.sum():,} px')
    im = axes[1, 0].imshow(prob_full, cmap='inferno', vmin=0, vmax=1)
    axes[1, 0].scatter(xs, ys, c='cyan', s=15, edgecolor='k')
    axes[1, 0].set_title(f'4. Outfall decay probability (scale {decay_scale:.0f} m)')
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


def main():
    qc_dir.mkdir(parents=True, exist_ok=True)
    for src_path in sorted(input_dir.glob('*_TURB.tif')):
        if src_path.name.startswith('._'):     # macOS metadata files on external drives
            continue
        mask_path = output_dir / src_path.name.replace('_TURB.tif', '_TURB_mask.tif')
        prob_path = output_dir / src_path.name.replace('_TURB.tif', '_TURB_prob.tif')
        qc_path = qc_dir / src_path.name.replace('_TURB.tif', '_TURB_qc.png')

        if mask_path.exists() and prob_path.exists() and qc_path.exists():
            print(f'Skipping {src_path.name} (already done)')
            continue

        print(f'Processing {src_path.name}')
        process_image(src_path, mask_path, prob_path, qc_path)

    print('Done.')


if __name__ == '__main__':
    main()
