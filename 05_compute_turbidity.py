"""
Compute Nechad turbidity (FNU) over water from Sentinel-2 reflectance mosaics.

Run after dn_to_reflectance.py.

    NDWI = (B03 - B08) / (B03 + B08)    -> water mask (NDWI >= 0.01) built ONCE from
                                           the clear reference date WATER_MASK_DATE and
                                           applied to every date; other pixels -> NaN
    T    = A * rho / (1 - rho / C)      turbidity in FNU, rho = B04 (665 nm)
                                        A = 610.94, C = 0.2324: Sentinel-2 665 nm coefficients,
                                        Nechad et al. (2016), ESA Living Planet Symposium, Table 3
                                        (B = 0). Pixels with rho < 0 or rho >= C
                                        (formula undefined / saturated, e.g. cloud) -> NaN.

Note: rho here is Sen2Cor L2A surface reflectance, used as an approximation of
water-leaving reflectance (rho_w). Sen2Cor is a land-oriented correction, so
glint/haze residuals feed straight into T. A and C are wavelength-specific.
Earlier versions used A = 247.10, C = 0.1697 (Nechad et al., 2009, at 657.5 nm,
the RapidEye red band); ACOLITE's Sentinel-2 B4 calibration is A = 366.14,
C = 0.19563.

Band order assumed in the mosaics: 3 = B03, 4 = B04, 8 = B08.

NDTI and NDCI are no longer computed here; the earlier version that also wrote
05_NDTI/ and 05_NDCI/ is in _archive/05_compute_ndti_ndci_full.py.

A fixed mask avoids hazy/glinty dates (e.g. 2022-05-13) where haze raises B03 and
B08 alike, NDWI collapses to ~0 and open water gets dropped as land. It does NOT
remove clouds or haze: pixels under cloud on other dates stay in and need a
separate cloud/haze screen. The mask is also saved as water_mask_<date>.tif.
"""

from pathlib import Path

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling

input_dir = Path("/home/jovyan/s2/03_mosaics_2022_2025_reflectance")
turb_dir = Path("/home/jovyan/s2/05_TURB")
turb_dir.mkdir(parents=True, exist_ok=True)

NDWI_THRESHOLD = 0.01
WATER_MASK_DATE = "20240825"   # clear, cloud-free reference scene

TURB_A = 610.94    # Nechad et al. (2016) Table 3, Sentinel-2 665 nm (B04)
TURB_C = 0.2324


def nechad_turbidity(rho):
    """Turbidity (FNU) from red reflectance; NaN where rho < 0 or rho >= C."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t = TURB_A * rho / (1.0 - rho / TURB_C)
    t[(rho < 0) | (rho >= TURB_C)] = np.nan
    return t


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
with rasterio.open(turb_dir.parent / f"water_mask_{WATER_MASK_DATE}.tif", "w", **mask_profile) as dst:
    dst.write(water.astype("uint8"), 1)

def water_on_grid(shape, transform, crs):
    """The water mask on a scene's own pixel grid. Mosaics normally share the
    reference grid; a scene with a different extent (e.g. a partial acquisition
    clipped to a smaller area) gets the mask resampled onto its grid
    (nearest neighbour; outside the reference extent = not water)."""
    if shape == water.shape and transform == mask_profile["transform"]:
        return water
    out = np.zeros(shape, dtype="uint8")
    reproject(water.astype("uint8"), out,
              src_transform=mask_profile["transform"], src_crs=mask_profile["crs"],
              dst_transform=transform, dst_crs=crs, resampling=Resampling.nearest)
    return out.astype(bool)


for src_path in sorted(input_dir.glob("*.tif")):
    turb_path = turb_dir / src_path.name.replace("_reflectance.tif", "_TURB.tif")

    if turb_path.exists():
        print(f"Skipping {src_path.name} (already done)")
        continue

    print(f"Processing {src_path.name}")

    with rasterio.open(src_path) as src:
        green = src.read(3)   # B03
        red = src.read(4)     # B04
        nir = src.read(8)     # B08
        profile = src.profile

    scene_water = water_on_grid(green.shape, profile["transform"], profile["crs"])
    if scene_water is not water:
        print(f"  note: grid {green.shape} differs from the water-mask grid {water.shape}; "
              f"mask resampled onto this scene")
    land = ~scene_water | ~np.isfinite(green) | ~np.isfinite(red)

    turb = nechad_turbidity(red)
    turb[land] = np.nan

    profile.update(count=1, dtype="float32", nodata=np.nan, compress="deflate")

    with rasterio.open(turb_path, "w", **profile) as dst:
        dst.write(turb.astype("float32"), 1)

print("Done.")
