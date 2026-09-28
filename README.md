# Sentinel-2 turbidity workflow (San Diego / Tijuana River coast)

<img width="200" height="400" alt="image" src="https://github.com/user-attachments/assets/641d18f3-4dce-44ad-b5d0-2a4f4a45f009" />

Downloads Sentinel-2 L2A scenes, mosaics and clips them to the study area,
converts them to surface reflectance, computes water-quality indices (NDTI, NDCI),
checks them with QC figures, and segments turbid plumes near outfalls.

Scripts run in order. Paths are set near the top of each script and assume the
`/home/jovyan/s2/` layout; edit them if your folders differ.

## Steps

| Step | Script | Input → output | What it does |
|---|---|---|---|
| 01 | `01_save.ipynb` | Copernicus Data Space → S3 bucket | Searches the CDSE STAC catalogue for L2A scenes over the AOI and date range, and copies the `.SAFE` products to S3 with rclone. Set the AOI, dates, remotes and `DST_PREFIX` in the first cell (`LIMIT` > 0 for a small test). |
| 02 | `02_merge_SAFEtoTIF.sh` | `.SAFE` folders → `tifs_2022_2025/` | Stacks 11 bands of each granule into one GeoTIFF (UInt16 DN, 20 m bands + 10 m bands). Run from the folder that contains the `.SAFE` products. |
| 03 | `03_merge_tifs_subset.py` | `tifs_2022_2025/` → `mosaics_2022_2025/` | Mosaics granules from the same acquisition time and clips them to a lon/lat polygon. Output: `mosaic_<YYYYMMDDTHHMMSS>_clipped.tif`. |
| 04 | `04_dn_to_reflectance.py` | mosaics → `03_mosaics_2022_2025_reflectance/` | DN → bottom-of-atmosphere reflectance: `(DN − 1000) / 10000`. NoData (DN 0) → NaN. |
| 05 | `05_compute_ndti_ndci.py` | reflectance → `05_NDTI/`, `05_NDCI/`, `water_mask_<date>.tif` | Computes NDTI and NDCI over water using one fixed water mask (see below). |
| 05b | `05b_qc_figures.py` | reflectance + NDTI → `06_QC/` | QC figures per date (true colour, spectral profiles at two points, NDTI map and histogram), an all-dates overview and `QC_stats.csv`. |
| 06 | `06_isegprob_ndti*.py`, `06_isegprob_ndti_tuning.ipynb` | NDTI → plume masks | Smooths and thresholds NDTI, cleans the mask with morphology, and keeps turbid clusters near outfalls (`shapefiles/Outflow.shp`). Tune on one scene in the notebook first. |

Step 03 example:

```bash
python 03_merge_tifs_subset.py \
  ../tifs_2022_2025 \
  ../mosaics_2022_2025 \
  --wkt "POLYGON ((-117.35 32.3, -117.02 32.3, -117.02 32.75, -117.35 32.75, -117.35 32.3))"
```

## Band order in the mosaics

Created in step 02 and assumed by every later step (1-based, as read by rasterio):

| 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
|---|---|---|---|---|---|---|---|---|---|---|
| B01 | B02 | B03 | B04 | B05 | B06 | B07 | B08 | B8A | B11 | B12 |

## Indices

- **NDWI** = (B03 − B08) / (B03 + B08). Used only to build the water mask (water = NDWI ≥ 0.01).
- **NDTI** = (B04 − B03) / (B04 + B03). Turbidity: sediment raises red, so higher = more turbid.
- **NDCI** = (B05 − B04) / (B05 + B04). Chlorophyll.

## Notes

- **Reflectance offset.** ESA added `BOA_ADD_OFFSET = −1000` with processing baseline
  04.00 (25 Jan 2022). The scenes downloaded in step 01 come from the reprocessed
  Collection-1 archive (baseline ≥ 05.00), which carries the offset for *all* dates,
  including earlier ones, so step 04 applies it everywhere. Sanity check: raw blue DN over
  dark water sits at about 1000–1200. Only use `OFFSET = 0` for genuine baseline 03.xx
  products.
- **Fixed water mask.** Step 05 builds the mask once from a clear reference date
  (`WATER_MASK_DATE`, currently `20240825`) and applies it to every date. Computing NDWI
  per date dropped open water on hazy scenes (2022-05-13 lost 41%). If you change the
  date, also update `WATER_MASK_LABEL` in `05b_qc_figures.py`.
- **Clouds and haze are not masked.** Pixels under cloud or haze still get NDTI values,
  biased toward 0. Screen or exclude these dates (e.g. median B08 over water > ~0.05)
  before comparing NDTI through time.
- **Very dark winter scenes.** Red reflectance over clear water can be ~0 or slightly
  negative after atmospheric correction, which pushes NDTI below −1. Treat those values as
  noise.
- **Tile seams.** Some dates show a straight jump in NDTI where two granules meet in
  the mosaic.
- **macOS external drives** create `._*.tif` files next to each GeoTIFF. The scripts glob
  `*.tif`, so delete those or filter them out before running.
