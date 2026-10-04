"""
s2_mask_overlays.py - S2 plume mask drawn over the true-colour image, one PNG per date

For every date that step 06 kept (status 'ok' in 06_summary.csv; hazy dates are
left out), writes a two-panel figure to OUT_DIR:
  left:  true colour (B4-B3-B2), 2nd-98th percentile stretch of that scene
  right: the same image with the plume mask from step 06 in red
Both panels mark the outfalls (cyan triangles, labelled with their names: TJ =
Tijuana River mouth, PB = Punta Bandera).
Reflectance comes from step 04 if present, otherwise from the step 03 DN mosaic
converted the same way as step 04: (DN - 1000) / 10000.
PNGs that already exist are skipped, so the script can be rerun.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
import geopandas as gpd
import matplotlib
matplotlib.use('Agg')                      # write PNGs without a display
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# --- Paths: TIDE cluster if /home/jovyan/s2 exists, otherwise the local Mac drive ---
TIDE = Path('/home/jovyan/s2')
if TIDE.exists():
    PLUME_DIR = TIDE / '06_TURB_plumes'                           # step 06 output
    REFLECTANCE_DIR = TIDE / '03_mosaics_2022_2025_reflectance'   # step 04 output
    DN_DIR = TIDE / '02_mosaics_2022_2025'                        # step 03 output
    OUTFALLS = TIDE / 'shapefiles' / 'Outflow.shp'
else:
    TJ = Path('/Volumes/External/TJ')
    PLUME_DIR = TJ / '010_optical' / '06_TURB_plumes'
    REFLECTANCE_DIR = TJ / '010_optical' / '03_mosaics_2022_2025_reflectance'   # not on the Mac
    DN_DIR = TJ / '010_optical' / '02_mosaics_2022_2025'
    OUTFALLS = TJ / '015_shapefiles' / 'outflow_PB_TJ' / 'Outflow.shp'
SUMMARY_CSV = PLUME_DIR / '06_summary.csv'
OUT_DIR = Path(__file__).resolve().parent / 's2_mask_overlays'


def read_true_colour(date):
    """B4-B3-B2 reflectance stack for one date, plus its grid; None if not found."""
    refl = REFLECTANCE_DIR / f'mosaic_{date}_reflectance.tif'
    dn = DN_DIR / f'mosaic_{date}_clipped.tif'
    if refl.exists():
        with rasterio.open(refl) as src:
            return np.dstack([src.read(b) for b in (4, 3, 2)]), src.transform, src.crs
    if dn.exists():
        with rasterio.open(dn) as src:
            bands = [src.read(b).astype(np.float32) for b in (4, 3, 2)]
            transform, crs = src.transform, src.crs
        return np.dstack([np.where(b == 0, np.nan, (b - 1000.0) / 10000.0) for b in bands]), transform, crs
    return None


def make_figure(row, outfalls, out_path):
    date = row['date']
    found = read_true_colour(date)
    if found is None:
        print(f'  no image for {date}')
        return
    rgb, transform, crs = found
    lo, hi = np.nanpercentile(rgb, [2, 98])
    rgb = np.nan_to_num(np.clip((rgb - lo) / (hi - lo), 0, 1), nan=1.0)

    with rasterio.open(PLUME_DIR / f'mosaic_{date}_TURB_mask.tif') as src:
        mask = src.read(1) == 1

    pts = outfalls.to_crs(crs)
    cols, rows = ~transform * (pts.geometry.x.values, pts.geometry.y.values)   # map -> pixel
    names = pts['featureNam'].tolist() if 'featureNam' in pts else [''] * len(pts)

    fig, axes = plt.subplots(1, 2, figsize=(13, 9.5), layout='constrained')
    for ax in axes:
        ax.imshow(rgb, interpolation='nearest')
        ax.scatter(cols, rows, marker='^', c='cyan', s=110, edgecolor='k', linewidth=1, zorder=4)
        for c, r, name in zip(cols, rows, names):
            ax.annotate(name, (c, r), xytext=(9, 4), textcoords='offset points', fontsize=11,
                        fontweight='bold', color='cyan', zorder=5,
                        path_effects=[pe.withStroke(linewidth=3, foreground='black')])
        ax.set_xticks([])
        ax.set_yticks([])
    axes[0].set_title('True colour (B4-B3-B2)')

    a = axes[1]
    a.imshow(np.ma.masked_where(~mask, mask), cmap=ListedColormap(['red']), alpha=0.45,
             interpolation='nearest')
    if mask.any():
        a.contour(mask, levels=[0.5], colors='red', linewidths=0.7)
    a.set_title(f'S2 plume mask: {row["plume_area_km2"]:.1f} km2')
    a.legend(handles=[Patch(color='red', alpha=0.6, label='S2 plume mask'),
                      Line2D([], [], marker='^', ls='', markerfacecolor='cyan', markeredgecolor='k',
                             markersize=9, label='outfall')], loc='lower left', fontsize=9)

    when = f'{date[:4]}-{date[4:6]}-{date[6:8]} {date[9:11]}:{date[11:13]} UTC'
    fig.suptitle(f'{when}   Sentinel-2', fontsize=13)
    fig.savefig(out_path, dpi=110, bbox_inches='tight')
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    outfalls = gpd.read_file(OUTFALLS)
    summary = pd.read_csv(SUMMARY_CSV, dtype={'date': str}).sort_values('date')
    summary = summary[summary['status'] == 'ok']           # hazy dates have no mask
    print(f'{len(summary)} dates with a plume mask')
    for _, row in summary.iterrows():
        out_path = OUT_DIR / f'S2_{row["date"]}_overlay.png'
        if out_path.exists():
            continue
        print(f'Plotting {row["date"]}')
        make_figure(row, outfalls, out_path)
    print('Done.')


if __name__ == '__main__':
    main()
