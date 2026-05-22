#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------
# Script para recorte y remuestreo VPD
# Vapor Pressure Deficit
# Interpolación recomendada: bilinear
# ------------------------------------

# Directorios
IN_DIR="/mnt/disco6tb/MHBA60/data/A1_RAW/VPD_v1.1"
temp_DIR="/mnt/disco6tb/MHBA60/data/A2_TEMP"
OUT_DIR="/mnt/disco6tb/MHBA60/data/A2_TEMP/VPD"

# Fichero de rejilla 0.25°
GRID_FILE="${temp_DIR}/grid_0.25.txt"

# Crear salida
mkdir -p "$OUT_DIR"

# Procesar todos los archivos TerraClimate VPD
for in_nc in "$IN_DIR"/TerraClimate_vpd_*.nc; do

  base=$(basename "$in_nc" .nc)
  out_nc="$OUT_DIR/${base}_0.25deg_bil.nc"

  echo "→ Procesando $base"
  echo "  Recorte global + remuestreo bilinear a 0.25°"

  cdo -O remapbil,"$GRID_FILE" \
      -sellonlatbox,-180,180,-90,90 \
      "$in_nc" "$out_nc"

done

cdo -O mergetime \
  /mnt/disco6tb/MHBA60/data/A2_TEMP/VPD/TerraClimate_vpd_{2003..2024}_0.25deg_bil.nc \
  /mnt/disco6tb/MHBA60/data/A2_TEMP/vpd_2003-2024_0.25deg_bil.nc



echo "Proceso finalizado."
