#!/usr/bin/env bash
set -euo pipefail

IN_DIR="/mnt/disco6tb/MHBA60/data/A1_RAW/MOD13C2_NDVI"
OUT_DIR="/mnt/disco6tb/MHBA60/data/A2_TEMP/MOD13C2_NDVI"

mkdir -p "$OUT_DIR"

shopt -s nullglob
files=("$IN_DIR"/MOD13C2.A*.hdf)

if [ ${#files[@]} -eq 0 ]; then
  echo "No se encontraron archivos .hdf en: $IN_DIR"
  exit 1
fi

for f in "${files[@]}"; do
  base=$(basename "$f" .hdf)
  out="$OUT_DIR/${base}.nc"

  echo "Procesando: $base"

  # Busca automáticamente el subdataset NDVI
  sds=$(gdalinfo "$f" | awk -F= '
    /SUBDATASET_[0-9]+_NAME=/ {
      name=$2
      low=tolower(name)
      if (low ~ /ndvi/ &&
          low !~ /evi/ &&
          low !~ /qa/ &&
          low !~ /quality/ &&
          low !~ /pixel reliability/) {
        print name
        exit
      }
    }'
  )

  if [ -z "$sds" ]; then
    echo "  -> No se encontró subdataset NDVI en $f"
    continue
  fi

  gdal_translate -of netCDF "$sds" "$out"
  echo "  -> Guardado en: $out"
done

echo "Conversión terminada."