"""
s1_s2_pairs_viewer.py - HTML page showing each matched S1/S2 pair side by side

Reads matched_pairs_with_areas.csv from step 07, keeps the pairs that have both an
S1 and an S2 plume mask and an S2 rating of 3 (the pairs used in
s1_s2_agreement_figure.py), makes any missing overlay PNGs with s1_mask_overlays.py and
s2_mask_overlays.py, and writes two pages next to this script:
  s1_s2_pairs.html        reads the full-size PNGs from s1_mask_overlays/ and
                          s2_mask_overlays/, so keep those folders next to it
  s1_s2_pairs_standalone.html  self-contained: every image is embedded in the file at
                          full resolution, so it can be shared on its own (e.g. on
                          Google Drive; viewers download it and open it in a browser)
Open either in a browser. The left/right arrow keys (or j / k) step through the pairs.
"""

import base64
import html
import importlib
import re
import sys
from pathlib import Path

import pandas as pd
import geopandas as gpd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
s1 = importlib.import_module('s1_mask_overlays')
s2 = importlib.import_module('s2_mask_overlays')

TIDE = Path('/home/jovyan/s2')
FIG07 = TIDE / '07_figures' if TIDE.exists() else Path('/Volumes/External/TJ/010_optical/07_figures')
PAIRS_CSV = FIG07 / 'matched_pairs_with_areas.csv'
OUT_HTML = HERE / 's1_s2_pairs.html'
OUT_HTML_STANDALONE = HERE / 's1_s2_pairs_standalone.html'


def load_pairs():
    """Matched pairs with both masks and S2 rating 3, with S1 scene and S2 date keys."""
    p = pd.read_csv(PAIRS_CSV)
    p = p[p['S1_File'].notna() & p['S2_File'].notna() & (p['S2_Rating'] == 3)].copy()
    p['s1_scene'] = p['S1_File'].map(lambda f: Path(f).name.split('_pre_')[0])
    p['s2_date'] = p['S2_File'].map(lambda f: re.search(r'\d{8}T\d{6}', Path(f).name).group(0))
    p['S1_Time'] = pd.to_datetime(p['S1_Time'], utc=True)
    p['S2_Time'] = pd.to_datetime(p['S2_Time'], utc=True)
    return p.sort_values('S1_Time').reset_index(drop=True)


def make_missing_overlays(pairs):
    """Create the overlay PNGs this page needs, if they are not there yet."""
    s1.OUT_DIR.mkdir(parents=True, exist_ok=True)
    s2.OUT_DIR.mkdir(parents=True, exist_ok=True)
    outfalls = gpd.read_file(s1.OUTFALLS)
    summary = pd.read_csv(s2.SUMMARY_CSV, dtype={'date': str}).set_index('date', drop=False)
    for _, r in pairs.iterrows():
        png1 = s1.OUT_DIR / f'{r.s1_scene}_overlay.png'
        if not png1.exists():
            print(f'S1 overlay {r.S1_Time:%Y-%m-%d}')
            mask = next(s1.MASK_DIR.glob(f'{r.s1_scene}_pre_*_mask.tif'))
            s1.make_figure(s1.PRE_DIR / f'{r.s1_scene}_pre.tif', mask, outfalls, png1)
        png2 = s2.OUT_DIR / f'S2_{r.s2_date}_overlay.png'
        if not png2.exists():
            print(f'S2 overlay {r.S2_Time:%Y-%m-%d}')
            s2.make_figure(summary.loc[r.s2_date], outfalls, png2)


def fmt(x, spec, unit=''):
    return '-' if pd.isna(x) else f'{x:{spec}}{unit}'


def embedded(rel_path):
    """The PNG, unchanged, as a data: URI for the standalone page."""
    return 'data:image/png;base64,' + base64.b64encode((HERE / rel_path).read_bytes()).decode('ascii')


def write_html(pairs, out_path, embed=False):
    rows, cards = [], []
    for i, r in pairs.iterrows():
        anchor = f'p{i}'
        gap_h = (r.S2_Time - r.S1_Time).total_seconds() / 3600
        stats = (f'Within the shared area: S1 {fmt(r.S1_Area_km2, ".1f", " km²")} · S2 {fmt(r.S2_Area_km2, ".1f", " km²")} · '
                 f'centroid distance {fmt(r.Centroid_Dist_m / 1000, ".1f", " km")} · '
                 f'overlap (IoU) {fmt(r.IoU_percent, ".1f", "%")} · S2 {gap_h:.1f} h after S1')
        rows.append(f'<tr><td><a href="#{anchor}">{r.S1_Time:%Y-%m-%d}</a></td>'
                    f'<td>{fmt(r.S1_Area_km2, ".1f")}</td><td>{fmt(r.S2_Area_km2, ".1f")}</td>'
                    f'<td>{fmt(r.Centroid_Dist_m / 1000, ".1f")}</td><td>{fmt(r.IoU_percent, ".1f")}</td></tr>')
        img1 = f's1_mask_overlays/{html.escape(r.s1_scene)}_overlay.png'
        img2 = f's2_mask_overlays/S2_{r.s2_date}_overlay.png'
        if embed:                      # browsers will not open data: images in a new tab
            pic1 = f'<img src="{embedded(img1)}" alt="S1 {r.S1_Time:%Y-%m-%d}">'
            pic2 = f'<img src="{embedded(img2)}" alt="S2 {r.S2_Time:%Y-%m-%d}">'
        else:
            pic1 = f'<a href="{img1}" target="_blank"><img src="{img1}" loading="lazy" alt="S1 {r.S1_Time:%Y-%m-%d}"></a>'
            pic2 = f'<a href="{img2}" target="_blank"><img src="{img2}" loading="lazy" alt="S2 {r.S2_Time:%Y-%m-%d}"></a>'
        cards.append(f'''
<section class="pair" id="{anchor}">
  <h2>{r.S1_Time:%Y-%m-%d} <span class="n">pair {i + 1} of {len(pairs)}</span></h2>
  <p class="stats">{stats}</p>
  <div class="imgs">
    <figure><figcaption>Sentinel-1 · {r.S1_Time:%H:%M} UTC</figcaption>
      {pic1}</figure>
    <figure><figcaption>Sentinel-2 · {r.S2_Time:%H:%M} UTC</figcaption>
      {pic2}</figure>
  </div>
  <p class="nav"><a href="#top">↑ list of pairs</a></p>
</section>''')

    page = f'''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>S1 / S2 plume pairs</title>
<style>
  body {{ font-family: -apple-system, Helvetica, Arial, sans-serif; margin: 0 auto; max-width: 1900px;
         padding: 16px; background: #fafafa; color: #222; }}
  h1 {{ margin: 0 0 4px; font-size: 22px; }}
  .lede {{ color: #555; margin: 0 0 16px; }}
  table {{ border-collapse: collapse; font-size: 14px; margin-bottom: 28px; }}
  th, td {{ padding: 4px 12px; text-align: right; border-bottom: 1px solid #ddd; }}
  th:first-child, td:first-child {{ text-align: left; }}
  th {{ background: #eee; position: sticky; top: 0; }}
  .pair {{ background: #fff; border: 1px solid #ddd; border-radius: 6px; padding: 12px 16px;
           margin-bottom: 24px; scroll-margin-top: 8px; }}
  .pair h2 {{ margin: 0; font-size: 19px; }}
  .n {{ font-weight: normal; color: #888; font-size: 14px; margin-left: 8px; }}
  .stats {{ margin: 4px 0 10px; color: #444; font-size: 14px; }}
  .imgs {{ display: flex; gap: 12px; flex-wrap: wrap; }}
  figure {{ flex: 1 1 600px; margin: 0; }}
  figcaption {{ font-weight: 600; margin-bottom: 4px; }}
  img {{ width: 100%; height: auto; border: 1px solid #ccc; }}
  .nav {{ margin: 8px 0 0; font-size: 13px; }}
</style></head>
<body id="top">
<h1>S1 / S2 plume pairs</h1>
<p class="lede">{len(pairs)} matched pairs with an S1 and an S2 plume mask (S2 dates rated 3).
Left: S1 mask on VV backscatter. Right: S2 mask on true colour.
{"" if embed else "Click an image to open it full size; "}arrow keys (or j / k) step through the pairs.<br>
Areas, centroid distance and overlap (IoU) in the table and above each pair are from step 07 and are measured
only within the area covered by both sensors (the S2 image extent). The S1 figure title gives the mask area
over the whole S1 scene, which extends further south and west, so it can be larger.</p>
<table>
<tr><th>Date</th><th>S1 area in shared area (km²)</th><th>S2 area (km²)</th><th>Centroid distance (km)</th><th>IoU (%)</th></tr>
{''.join(rows)}
</table>
{''.join(cards)}
<script>
  const pairs = [...document.querySelectorAll('.pair')];
  function current() {{
    let best = 0;
    pairs.forEach((p, i) => {{ if (p.getBoundingClientRect().top <= 60) best = i; }});
    return best;
  }}
  document.addEventListener('keydown', e => {{
    if (e.target.tagName === 'INPUT') return;
    let i = current();
    if (e.key === 'ArrowRight' || e.key === 'j') i = Math.min(i + 1, pairs.length - 1);
    else if (e.key === 'ArrowLeft' || e.key === 'k') i = Math.max(i - 1, 0);
    else return;
    e.preventDefault();
    pairs[i].scrollIntoView({{ behavior: 'smooth' }});
  }});
</script>
</body></html>
'''
    out_path.write_text(page, encoding='utf-8')
    print(f'Wrote {out_path} ({out_path.stat().st_size / 1e6:.1f} MB)')


def main():
    pairs = load_pairs()
    print(f'{len(pairs)} matched pairs with both masks (S2 rating 3)')
    make_missing_overlays(pairs)
    write_html(pairs, OUT_HTML)
    write_html(pairs, OUT_HTML_STANDALONE, embed=True)


if __name__ == '__main__':
    main()
