find . -path "*/GRANULE/*/IMG_DATA/R20m/*_B01_20m.jp2" -print0 | \
xargs -0 -n 1 -P 4 bash -c '
  i="$1"

  safe_dir=$(echo "$i" | sed -n "s#.*/\([^/]*\.SAFE\)/GRANULE/.*#\1#p")
  product="${safe_dir%.SAFE}"
  granule=$(echo "$i" | sed -n "s#.*/GRANULE/\([^/]*\)/IMG_DATA/.*#\1#p")

  outdir="$HOME/s2/tifs_2024"
  mkdir -p "$outdir"

  outfile="$outdir/${product}_${granule}_L2A_e.tif"
  vrt="${outfile}.vrt"

  B01="$i"
  B05="${i/_B01_20m.jp2/_B05_20m.jp2}"
  B06="${i/_B01_20m.jp2/_B06_20m.jp2}"
  B07="${i/_B01_20m.jp2/_B07_20m.jp2}"
  B8A="${i/_B01_20m.jp2/_B8A_20m.jp2}"
  B11="${i/_B01_20m.jp2/_B11_20m.jp2}"
  B12="${i/_B01_20m.jp2/_B12_20m.jp2}"

  B02="${i/IMG_DATA\/R20m/IMG_DATA\/R10m}"
  B02="${B02/_B01_20m.jp2/_B02_10m.jp2}"

  B03="${i/IMG_DATA\/R20m/IMG_DATA\/R10m}"
  B03="${B03/_B01_20m.jp2/_B03_10m.jp2}"

  B04="${i/IMG_DATA\/R20m/IMG_DATA\/R10m}"
  B04="${B04/_B01_20m.jp2/_B04_10m.jp2}"

  B08="${i/IMG_DATA\/R20m/IMG_DATA\/R10m}"
  B08="${B08/_B01_20m.jp2/_B08_10m.jp2}"

  gdalbuildvrt -separate -q "$vrt" \
    "$B01" "$B02" "$B03" "$B04" "$B05" "$B06" "$B07" "$B08" "$B8A" "$B11" "$B12" && \
  gdal_translate -q -of GTiff -ot UInt16 \
    -co TILED=YES -co COMPRESS=DEFLATE -co BIGTIFF=YES \
    "$vrt" "$outfile"

  rm -f "$vrt"
' _