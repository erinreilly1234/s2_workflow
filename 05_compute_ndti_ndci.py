"""
Compute NDTI and NDCI over water from Sentinel-2 reflectance mosaics.

Run after dn_to_reflectance.py.

    NDWI = (B03 - B08) / (B03 + B08)    -> water mask (NDWI >= 0.01) built ONCE from
                                           the clear reference date WATER_MASK_DATE and
                                           applied to every date; other pixels -> NaN
    NDTI = (B04 - B03) / (B04 + B03)    turbidity
    NDCI = (B05 - B04) / (B05 + B04)    chlorophyll

Band order assumed in the mosaics: 2 = B02, 3 = B03, 4 = B04, 5 = B05, 8 = B08.

A fixed mask avoids hazy/glinty dates (e.g. 2022-05-13) where haze raises B03 and
B08 alike, NDWI collapses to ~0 and open water gets dropped as land. It does NOT
remove clouds or haze: pixels under cloud on other dates stay in and need a
separate cloud/haze screen. The mask is also saved as water_mask_<date>.tif.
"""

from pathlib import Path

import numpy as np
import rasterio

input_dir = Path("/home/jovyan/s2/03_mosaics_2022_2025_reflectance")
ndti_dir = Path("/home/jovyan/s2/05_NDTI")
ndci_dir = Path("/home/jovyan/s2/05_NDCI")
ndti_dir.mkdir(parents=True, exist_ok=True)
ndci_dir.mkdir(parents=True, exist_ok=True)

NDWI_THRESHOLD = 0.01
WATER_MASK_DATE = "20240825"   # clear, cloud-free reference scene


def normalized_difference(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (a - b) / (a + b)


ref_path = next(input_dir.glob(f"mosaic_{WATER_MASK_DATE}*_reflectance.tif"))
with rasterio.open(ref_path) as src:
    ndwi_ref = normalized_difference(src.read(3), src.read(8))
    mask_profile = src.profile
water = ndwi_ref >= NDWI_THRESHOLD   # NaN compares False -> not water
print(f"Water mask from {ref_path.name}: {water.sum():,} water pixels")

mask_profile.update(count=1, dtype="uint8", nodata=255, compress="deflate")
with rasterio.open(ndti_dir.parent / f"water_mask_{WATER_MASK_DATE}.tif", "w", **mask_profile) as dst:
    dst.write(water.astype("uint8"), 1)

for src_path in sorted(input_dir.glob("*.tif")):
    ndti_path = ndti_dir / src_path.name.replace("_reflectance.tif", "_NDTI.tif")
    ndci_path = ndci_dir / src_path.name.replace("_reflectance.tif", "_NDCI.tif")

    if ndti_path.exists() and ndci_path.exists():
        print(f"Skipping {src_path.name} (already done)")
        continue

    print(f"Processing {src_path.name}")

    with rasterio.open(src_path) as src:
        green = src.read(3)   # B03
        red = src.read(4)     # B04
        b05 = src.read(5)     # B05
        nir = src.read(8)     # B08
        profile = src.profile

    if green.shape != water.shape:
        raise ValueError(f"{src_path.name} grid {green.shape} != water mask grid {water.shape}")
    land = ~water | ~np.isfinite(green) | ~np.isfinite(red)

    ndti = normalized_difference(red, green)
    ndci = normalized_difference(b05, red)
    ndti[land] = np.nan
    ndci[land] = np.nan

    profile.update(count=1, dtype="float32", nodata=np.nan, compress="deflate")

    with rasterio.open(ndti_path, "w", **profile) as dst:
        dst.write(ndti.astype("float32"), 1)

    with rasterio.open(ndci_path, "w", **profile) as dst:
        dst.write(ndci.astype("float32"), 1)

print("Done.")
