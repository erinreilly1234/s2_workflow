#!/usr/bin/env python3
from __future__ import annotations

import re
from pathlib import Path
from collections import defaultdict

import rasterio
from rasterio.merge import merge
from rasterio.mask import mask

from shapely import wkt
from shapely.geometry import mapping
from shapely.ops import transform as shp_transform
from pyproj import Transformer


DT_RE = re.compile(r"MSIL2A_(\d{8}T\d{6})")


def extract_datetime_from_name(name: str) -> str | None:
    m = DT_RE.search(name)
    return m.group(1) if m else None


def clip_array_to_wkt(mosaic, transform, crs, polygon_wkt: str):
    """
    Clip a mosaic array to polygon_wkt (assumed EPSG:4326 lon/lat),
    reprojecting geometry to raster CRS if needed.
    Returns (clipped_array, clipped_transform).
    """
    geom_ll = wkt.loads(polygon_wkt)  # lon/lat polygon

    # Reproject polygon to raster CRS (if raster CRS differs from EPSG:4326)
    if crs is None:
        raise ValueError("Raster CRS is missing; cannot reproject study area polygon safely.")

    raster_crs = crs
    transformer = Transformer.from_crs("EPSG:4326", raster_crs, always_xy=True)
    geom_raster = shp_transform(transformer.transform, geom_ll)

    # rasterio.mask.mask expects shapes in raster CRS
    # But mask() operates on an open dataset, not raw arrays.
    # We'll write to an in-memory dataset using rasterio.io.MemoryFile.
    from rasterio.io import MemoryFile

    # Build a temporary dataset in memory
    meta = {
        "driver": "GTiff",
        "height": mosaic.shape[1],
        "width": mosaic.shape[2],
        "count": mosaic.shape[0],
        "dtype": mosaic.dtype,
        "crs": raster_crs,
        "transform": transform,
    }

    with MemoryFile() as memfile:
        with memfile.open(**meta) as ds:
            ds.write(mosaic)
            out_image, out_transform = mask(
                ds,
                [mapping(geom_raster)],
                crop=True,
                all_touched=False,  # set True if you want more inclusive edge pixels
                filled=True,
            )
            out_meta = ds.meta.copy()

    return out_image, out_transform, out_meta


def mosaic_and_clip_group(tif_paths: list[Path], out_path: Path, polygon_wkt: str) -> None:
    srcs = [rasterio.open(p) for p in tif_paths]
    try:
        mosaic, out_transform = merge(srcs)

        # template metadata from first
        out_meta = srcs[0].meta.copy()
        out_meta.update(
            {
                "height": mosaic.shape[1],
                "width": mosaic.shape[2],
                "transform": out_transform,
                "count": mosaic.shape[0],
                "compress": "deflate",
                "tiled": True,
                "BIGTIFF": "IF_SAFER",
            }
        )

        # Clip to study area
        clipped, clipped_transform, tmp_meta = clip_array_to_wkt(
            mosaic, out_transform, out_meta.get("crs"), polygon_wkt
        )

        out_meta.update(
            {
                "height": clipped.shape[1],
                "width": clipped.shape[2],
                "transform": clipped_transform,
                "count": clipped.shape[0],
                "dtype": clipped.dtype,
            }
        )

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with rasterio.open(out_path, "w", **out_meta) as dest:
            dest.write(clipped)

    finally:
        for s in srcs:
            s.close()


def main(in_dir: str, out_dir: str, polygon_wkt: str, pattern: str = "*.tif") -> None:
    in_path = Path(in_dir)
    out_path = Path(out_dir)

    tifs = sorted(in_path.glob(pattern))
    if not tifs:
        raise SystemExit(f"No files matching {pattern} in {in_path}")

    groups: dict[str, list[Path]] = defaultdict(list)
    skipped = []

    for tif in tifs:
        dt = extract_datetime_from_name(tif.name)
        if dt is None:
            skipped.append(tif.name)
            continue
        groups[dt].append(tif)

    if skipped:
        print("Skipped (no datetime found):")
        for s in skipped:
            print("  ", s)

    print(f"Found {len(groups)} datetime groups")

    for dt, paths in sorted(groups.items()):
        out_file = out_path / f"mosaic_{dt}_clipped.tif"
        print(f"\nMerging {len(paths)} tiles for {dt} -> clip -> {out_file}")
        mosaic_and_clip_group(paths, out_file, polygon_wkt)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Group GeoTIFFs by datetime, mosaic each group, and clip to a WKT polygon.")
    ap.add_argument("in_dir", help="Folder containing .tif tiles")
    ap.add_argument("out_dir", help="Folder to write mosaicked + clipped outputs")
    ap.add_argument(
        "--wkt",
        required=True,
        help="Study area polygon in WKT (assumed lon/lat EPSG:4326)",
    )
    ap.add_argument("--pattern", default="*.tif", help="Glob pattern (default: *.tif)")
    args = ap.parse_args()

    main(args.in_dir, args.out_dir, args.wkt, args.pattern)