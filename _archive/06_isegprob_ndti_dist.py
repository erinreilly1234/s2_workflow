"""
06_isegprob_ndti.py

Run after 05_compute_ndti_ndci.py.

Adapts the SAR isegprob pipeline (threshold -> morphological cleanup ->
outfall-distance cluster selection) to the Sentinel-2 NDTI rasters written by
step 05:
  - low-pass filter with valid-data weighting (land pixels are already NaN
    from the NDWI mask in step 05)
  - NDTI >= threshold marks turbid water
  - morphological open/close + small-object/hole removal + dilation
  - rasterize outfall points, Euclidean distance transform (in meters)
  - connected-component labeling; a cluster is kept if any part of it lies
    within max_distance_m of an outfall

Distances are computed from the raster's own pixel size, so this works for the
native UTM mosaics produced by step 03 (and for geographic rasters too).
Smoothing window, morphology radii and minimum areas are set in meters and
converted to pixels per image, so the same settings work for SAR and optical
rasters with different pixel sizes.

Outputs per scene: a selected-cluster mask (_NDTI_mask.tif), a map of each
cluster's closest distance to an outfall in meters (_NDTI_dist.tif), and a
6-panel QC figure (raw NDTI, threshold, morphology, outfall distance, clusters,
final mask) in output_dir/qc_pngs.
"""

import math
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter, distance_transform_edt, label, minimum
from skimage.morphology import (
    binary_opening, binary_closing, binary_dilation,
    remove_small_objects, remove_small_holes, disk
)
import matplotlib
matplotlib.use('Agg')  # write PNGs without needing a display
import matplotlib.pyplot as plt
import geopandas as gpd
import rasterio
from rasterio.features import rasterize

# --- Parameters ---
input_dir = Path('/home/jovyan/s2/04_NDTI')             # NDTI output of step 05
output_dir = Path('/home/jovyan/s2/05_NDTI_isegprob')
qc_dir = output_dir / 'qc_pngs'                       # per-scene QC figures
make_qc_png = True                                    # set False to skip the figures
outfall_shapefile = '/home/jovyan/s2/shapefiles/Outflow.shp'

threshold = 0.07              # NDTI at or above this counts as turbid water
max_distance_m = 1300.0        # keep clusters that come within this distance of an outfall

# Spatial settings in meters / square meters, converted to pixels per image so
# the same values work for SAR and optical rasters with different pixel sizes.
# Defaults equal the original pixel values at a 10 m pixel.
smooth_window_m = 50.0         # low-pass filter window width
open_radius_m = 30.0           # opening: removes specks/thin features
close_radius_m = 50.0          # closing: bridges small gaps
min_hole_area_m2 = 80000.0     # fill holes smaller than this (3 ha)
min_object_area_m2 = 30000.0   # drop patches smaller than this (1 ha)
dilate_radius_m = 30.0         # final dilation buffer


def disk_px(radius_m, px_m):
    """Structuring element whose radius is radius_m meters (at least 1 pixel)."""
    return disk(max(1, round(radius_m / px_m)))


def area_px(area_m2, px_x, px_y):
    """Convert an area in square meters to a pixel count (at least 1 pixel)."""
    return max(1, round(area_m2 / (px_x * px_y)))


def save_qc_png(out_png, title, img, signal, mask, dist_m, outfall_r,
                labels, n_clusters, n_accepted, mask_selected):
    """Save a 6-panel QC figure for one scene."""
    clims = np.nanpercentile(img, [2, 98]) if np.any(~np.isnan(img)) else (0, 1)
    ys, xs = np.where(outfall_r == 1)
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    axes[0, 0].imshow(img, cmap='viridis', vmin=clims[0], vmax=clims[1])
    axes[0, 0].set_title('1. Raw NDTI')

    axes[0, 1].imshow(signal, cmap='gray')
    axes[0, 1].set_title(f'2. NDTI >= {threshold}\n{signal.sum():,} px')

    axes[0, 2].imshow(mask, cmap='gray')
    axes[0, 2].set_title(f'3. After morphology\n{mask.sum():,} px')

    im = axes[1, 0].imshow(dist_m / 1000.0, cmap='magma_r')
    axes[1, 0].contour(dist_m, levels=[max_distance_m], colors='cyan', linewidths=1.2)
    axes[1, 0].scatter(xs, ys, c='cyan', s=15, edgecolor='k')
    axes[1, 0].set_title(f'4. Distance to outfall (cyan = {max_distance_m:.0f} m cutoff)')
    fig.colorbar(im, ax=axes[1, 0], fraction=0.046, label='km')

    axes[1, 1].imshow(labels, cmap='nipy_spectral', interpolation='nearest')
    axes[1, 1].set_title(f'5. Clusters (n={n_clusters}, {n_accepted} accepted)')

    axes[1, 2].imshow(img, cmap='gray', vmin=clims[0], vmax=clims[1])
    overlay = np.ma.masked_where(~mask_selected, mask_selected)
    axes[1, 2].imshow(overlay, cmap='autumn', alpha=0.6, interpolation='nearest')
    axes[1, 2].contour(dist_m, levels=[max_distance_m], colors='cyan', linewidths=1.2)
    axes[1, 2].scatter(xs, ys, c='cyan', s=15, edgecolor='k')
    axes[1, 2].set_title(f'6. Final plume mask\n{mask_selected.sum():,} px')

    for ax in axes.ravel():
        ax.axis('off')
    fig.suptitle(title, fontsize=13)
    fig.tight_layout()
    fig.savefig(out_png, dpi=130, bbox_inches='tight')
    plt.close(fig)


def process_image(image_path, mask_path, dist_path, qc_path=None):
    # 1. Read NDTI raster
    with rasterio.open(image_path) as src:
        img = src.read(1).astype(np.float32)
        profile = src.profile
        transform = src.transform
        crs = src.crs
        h, w = src.height, src.width
        bounds = src.bounds

    # Pixel size in meters (UTM mosaics are already in meters)
    px_x, px_y = abs(transform.a), abs(transform.e)
    if crs.is_geographic:
        mid_lat = (bounds.bottom + bounds.top) / 2.0
        px_x *= 111320 * math.cos(math.radians(mid_lat))
        px_y *= 110540
    px_mean = (px_x + px_y) / 2.0
    win = max(1, round(smooth_window_m / px_mean))
    win += (win % 2 == 0)  # keep the window odd so it stays centered

    # 2. Low-pass filter with valid-data weighting
    valid = (~np.isnan(img)).astype(np.float32)
    filled = np.nan_to_num(img, nan=0.0)
    num = uniform_filter(filled, size=win, mode='nearest')
    den = uniform_filter(valid, size=win, mode='nearest')
    smooth = np.full_like(img, np.nan, dtype=np.float32)
    has_data = den > 0
    smooth[has_data] = num[has_data] / den[has_data]

    # 3. Threshold & morphological cleanup
    signal = (smooth >= threshold) & ~np.isnan(smooth)
    clean = binary_opening(signal, disk_px(open_radius_m, px_mean))
    clean = binary_closing(clean, disk_px(close_radius_m, px_mean))
    clean = remove_small_holes(clean, area_threshold=area_px(min_hole_area_m2, px_x, px_y))
    clean = remove_small_objects(clean, min_size=area_px(min_object_area_m2, px_x, px_y))
    mask = binary_dilation(clean, disk_px(dilate_radius_m, px_mean))

    # 4. Rasterize outfall points and compute distance (meters)
    outfalls = gpd.read_file(outfall_shapefile).to_crs(crs)
    outfall_r = rasterize(
        [(geom, 1) for geom in outfalls.geometry],
        out_shape=(h, w), transform=transform,
        fill=0, dtype='uint8'
    )
    if not outfall_r.any():
        print(f'  WARNING: no outfall falls inside {image_path.name}; skipping.')
        return
    dist_m = distance_transform_edt(outfall_r == 0, sampling=(px_y, px_x))

    dist_m = dist_m.astype(np.float32)

    # 5. Label clusters; keep those that come within max_distance_m of an outfall
    labels, n_clusters = label(mask)
    dist_map = np.full(mask.shape, np.nan, dtype=np.float32)  # NaN outside clusters
    mask_selected = np.zeros(mask.shape, dtype=bool)
    n_accepted = 0
    if n_clusters > 0:
        ids = np.arange(1, n_clusters + 1)
        min_d = np.asarray(minimum(dist_m, labels, index=ids), dtype=np.float32)
        lut = np.full(n_clusters + 1, np.nan, dtype=np.float32)
        lut[1:] = min_d
        dist_map = lut[labels]
        accepted = ids[min_d <= max_distance_m]
        n_accepted = len(accepted)
        mask_selected = np.isin(labels, accepted)

    # 6. Write selected-cluster mask
    mask_profile = profile.copy()
    mask_profile.update(driver='GTiff', dtype='uint8', count=1, nodata=0,
                        compress='deflate')
    with rasterio.open(mask_path, 'w', **mask_profile) as dst:
        dst.write(mask_selected.astype('uint8'), 1)

    # 7. Write each cluster's closest distance to an outfall (meters)
    dist_profile = profile.copy()
    dist_profile.update(driver='GTiff', dtype='float32', count=1, nodata=np.nan,
                        compress='deflate')
    with rasterio.open(dist_path, 'w', **dist_profile) as dst:
        dst.write(dist_map, 1)

    print(f'  Saved mask: {mask_path.name}')
    print(f'  Saved distance map: {dist_path.name}')

    # 8. QC figure
    if qc_path is not None:
        save_qc_png(qc_path, image_path.name, img, signal, mask, dist_m,
                    outfall_r, labels, n_clusters, n_accepted, mask_selected)
        print(f'  Saved QC png: {qc_path.name}')


def main():
    output_dir.mkdir(parents=True, exist_ok=True)
    if make_qc_png:
        qc_dir.mkdir(parents=True, exist_ok=True)
    for src_path in sorted(input_dir.glob('*_NDTI.tif')):
        mask_path = output_dir / src_path.name.replace('_NDTI.tif', '_NDTI_mask.tif')
        dist_path = output_dir / src_path.name.replace('_NDTI.tif', '_NDTI_dist.tif')
        qc_path = (qc_dir / src_path.name.replace('_NDTI.tif', '_NDTI_qc.png')
                   if make_qc_png else None)

        if mask_path.exists() and dist_path.exists() and (qc_path is None or qc_path.exists()):
            print(f'Skipping {src_path.name} (already done)')
            continue

        print(f'Processing {src_path.name}')
        process_image(src_path, mask_path, dist_path, qc_path)

    print('Done.')


if __name__ == '__main__':
    main()
