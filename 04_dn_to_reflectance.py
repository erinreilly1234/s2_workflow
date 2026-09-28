"""
Convert Sentinel-2 L2A mosaics from DN to BOA reflectance.

    reflectance = (DN + BOA_ADD_OFFSET) / BOA_QUANTIFICATION_VALUE

Values from MTD_MSIL2A.xml (processing baseline 05.10).

Note on the January 2022 change: ESA introduced BOA_ADD_OFFSET = -1000 with
processing baseline 04.00 (operational from 25 Jan 2022). Products still on
baseline 03.xx have no offset. The scenes used here come from the reprocessed
Collection-1 archive (baseline >= 05.00), which applies the offset to every
date, including acquisitions before 25 Jan 2022, so the offset is applied to
all files. Check: DN offset data has a floor near 1000 over dark water; if
you mix in any baseline 03.xx products, set OFFSET = 0.0 for those.

Output keeps the input filename; NoData pixels (DN 0) become NaN.
"""

from pathlib import Path

import numpy as np
import rasterio

input_dir = Path("/home/jovyan/s2/02_mosaics_2022_2025")
output_dir = Path("/home/jovyan/s2/03_mosaics_2022_2025_reflectance")
output_dir.mkdir(parents=True, exist_ok=True)

QUANTIFICATION = 10000.0
OFFSET = -1000.0

for src_path in sorted(input_dir.glob("*.tif")):
    
    output_name = src_path.stem.removesuffix("_clipped") + "_reflectance" + src_path.suffix
    out_path = output_dir / output_name

    if out_path.exists():
        print(f"Skipping {out_path.name} (already converted)")
        continue

    print(f"Converting {src_path.name}")

    with rasterio.open(src_path) as src:
        dn = src.read().astype("float32")
        profile = src.profile

    reflectance = (dn + OFFSET) / QUANTIFICATION
    reflectance[dn == 0] = np.nan

    profile.update(dtype="float32", nodata=np.nan, compress="deflate")

    with rasterio.open(out_path, "w", **profile) as dst:
        dst.write(reflectance)

print("Done.")
