"""
Compute NDTI and NDCI over water from Sentinel-2 reflectance mosaics.

Run after dn_to_reflectance.py.

    NDWI = (B03 - B08) / (B03 + B08)    -> pixels with NDWI < 0.1 are land (set to NaN)
    NDTI = (B04 - B03) / (B04 + B03)    turbidity
    NDCI = (B05 - B04) / (B05 + B04)    chlorophyll

Band order assumed in the mosaics: 2 = B02, 3 = B03, 4 = B04, 5 = B05, 8 = B08.
"""

from pathlib import Path

import numpy as np
import rasterio

input_dir = Path("/home/jovyan/s2/mosaics_2022_2025_reflectance")
ndti_dir = Path("/home/jovyan/s2/06_NDTI")
ndci_dir = Path("/home/jovyan/s2/06_NDCI")
ndti_dir.mkdir(parents=True, exist_ok=True)
ndci_dir.mkdir(parents=True, exist_ok=True)

NDWI_THRESHOLD = 0.01


def normalized_difference(a, b):
    with np.errstate(divide="ignore", invalid="ignore"):
        return (a - b) / (a + b)


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

    ndwi = normalized_difference(green, nir)
    land = ~(ndwi >= NDWI_THRESHOLD)   # NDWI < 0.1, plus NaN/NoData pixels

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
