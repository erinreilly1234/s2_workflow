"""
s1_mask_overlays.py - S1 plume mask drawn over the preprocessed S1 image, one PNG per date

For every preprocessed Sentinel-1 scene in PRE_DIR that has a plume mask in MASK_DIR,
writes a two-panel figure to OUT_DIR:
  left:  VV backscatter (sigma0, band 1 of the preprocessed GeoTIFF) in dB,
         2nd-98th percentile stretch of that scene; a 3x3 mean is applied for
         display only, to reduce speckle (the mask is drawn exactly as saved)
  right: the same image with the plume mask in red
Both panels mark the outfalls (cyan triangles, labelled with their names: TJ =
Tijuana River mouth, PB = Punta Bandera).
Files are matched by scene name: <scene>_pre.tif  <->  <scene>_pre_*_mask.tif.
PNGs that already exist are skipped, so the script can be rerun after adding scenes.
"""

import re
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter
import rasterio
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')                      # write PNGs without a display
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch

# --- Paths ---
TJ = Path('/Volumes/External/TJ')
PRE_DIR = TJ / '020_preprocessing'                      # <scene>_pre.tif (band 1 = VV sigma0)
MASK_DIR = TJ / '040_segoutput'                         # <scene>_pre_..._mask.tif (1 = plume)
OUTFALLS = TJ / '015_shapefiles' / 'outflow_PB_TJ' / 'Outflow.shp'
OUT_DIR = Path(__file__).resolve().parent / 's1_mask_overlays'


def scene_time(name):
    """'2022-01-01 01:50 UTC' from a Sentinel-1 file name."""
    t = re.search(r'_(\d{8}T\d{6})_', name).group(1)
    return f'{t[:4]}-{t[4:6]}-{t[6:8]} {t[9:11]}:{t[11:13]} UTC'


def area_km2(mask, transform, crs):
    """Area of the mask in km2, for a projected (m) or a lat/lon grid."""
    if crs.is_geographic:                  # pixel area shrinks with latitude
        lat = transform.f + transform.e * (np.arange(mask.shape[0]) + 0.5)
        row_area = abs(transform.a) * 111320 * np.cos(np.radians(lat)) * abs(transform.e) * 110540
        return float((mask.sum(axis=1) * row_area).sum() / 1e6)
    return float(mask.sum() * abs(transform.a * transform.e) / 1e6)


def make_figure(pre_path, mask_path, outfalls, out_path):
    with rasterio.open(pre_path) as src:
        sigma0 = src.read(1).astype(np.float32)
        transform, crs = src.transform, src.crs
    with rasterio.open(mask_path) as src:
        mask = src.read(1) == 1
    if mask.shape != sigma0.shape:
        raise ValueError(f'{mask_path.name} and {pre_path.name} are on different grids')

    valid = np.isfinite(sigma0) & (sigma0 > 0)
    total = uniform_filter(np.where(valid, sigma0, 0.0), size=3)
    count = uniform_filter(valid.astype(np.float32), size=3)
    with np.errstate(divide='ignore', invalid='ignore'):
        smooth = np.where(valid & (count > 0), total / count, np.nan)   # 3x3 mean, display only
        db = 10 * np.log10(smooth)
    lo, hi = np.nanpercentile(db, [2, 98])

    pts = outfalls.to_crs(crs)
    cols, rows = ~transform * (pts.geometry.x.values, pts.geometry.y.values)   # map -> pixel
    names = pts['featureNam'].tolist() if 'featureNam' in pts else [''] * len(pts)

    def mark_outfalls(ax):
        ax.scatter(cols, rows, marker='^', c='cyan', s=110, edgecolor='k', linewidth=1, zorder=4)
        for c, r, name in zip(cols, rows, names):
            ax.annotate(name, (c, r), xytext=(9, 4), textcoords='offset points', fontsize=11,
                        fontweight='bold', color='cyan', zorder=5,
                        path_effects=[pe.withStroke(linewidth=3, foreground='black')])

    fig, axes = plt.subplots(1, 2, figsize=(14, 7.5), layout='constrained')
    for ax in axes:
        im = ax.imshow(db, cmap='gray', vmin=lo, vmax=hi, interpolation='nearest')
        mark_outfalls(ax)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.colorbar(im, ax=list(axes), shrink=0.8, pad=0.01, label='VV backscatter (dB)')
    axes[0].set_title('VV backscatter')

    a = axes[1]
    a.imshow(np.ma.masked_where(~mask, mask), cmap=ListedColormap(['red']), alpha=0.45,
             interpolation='nearest')
    if mask.any():
        a.contour(mask, levels=[0.5], colors='red', linewidths=0.7)
    a.legend(handles=[Patch(color='red', alpha=0.6, label='S1 plume mask'),
                      Line2D([], [], marker='^', ls='', markerfacecolor='cyan', markeredgecolor='k',
                             markersize=9, label='outfall')], loc='lower left', fontsize=9)
    a.set_title(f'S1 plume mask: {area_km2(mask, transform, crs):.1f} km2')

    fig.suptitle(f'{scene_time(pre_path.name)}   {pre_path.name[:3]}', fontsize=13)
    fig.savefig(out_path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    outfalls = gpd.read_file(OUTFALLS)
    masks = {p.name.split('_pre_')[0]: p for p in MASK_DIR.glob('*_mask.tif') if not p.name.startswith('._')}
    pres = sorted(p for p in PRE_DIR.glob('*_pre.tif') if not p.name.startswith('._'))
    print(f'{len(pres)} preprocessed scenes, {len(masks)} masks')

    for pre_path in pres:
        scene = pre_path.name[:-len('_pre.tif')]
        if scene not in masks:
            print(f'No mask for {scene}')
            continue
        out_path = OUT_DIR / f'{scene}_overlay.png'
        if out_path.exists():
            continue
        print(f'Plotting {scene_time(scene + "_")}')
        make_figure(pre_path, masks[scene], outfalls, out_path)

    print('Done.')


if __name__ == '__main__':
    main()
