#!/usr/bin/env bash
set -euo pipefail

# ------------------------------------
# Preprocesado SMs GLEAM v4.2b
# Unión temporal + remuestreo a 0.25°
# ------------------------------------

# Directorios
IN_DIR="/mnt/disco6tb/MHBA60/data/A1_RAW/SMs_v4.2b"
OUT_DIR="/mnt/disco6tb/MHBA60/data/A2_TEMP"

# Rejilla destino 0.25°
GRID_025="${OUT_DIR}/grid_0.25.txt"

# Rejilla original GLEAM 0.1°
SRC_GRID="${OUT_DIR}/srcgrid_GLEAM_0.1deg.txt"

# Archivos de salida
FILELIST="${OUT_DIR}/filelist_SMs_GLEAM_v4.2b.txt"
MERGED_NC="${OUT_DIR}/SMs_2003_2024_GLEAM_v4.2b_MO.nc"
OUT_NC="${OUT_DIR}/SMs_2003_2024_GLEAM_v4.2b_MO_025deg_bil.nc"

mkdir -p "$OUT_DIR"

# Comprobar rejilla destino
if [[ ! -f "$GRID_025" ]]; then
  echo "ERROR: no existe la rejilla destino:"
  echo "$GRID_025"
  exit 1
fi

# Crear rejilla original GLEAM 0.1°
cat > "$SRC_GRID" << EOF
gridtype     = lonlat
xsize        = 3600
ysize        = 1800
xfirst       = -179.95
yfirst       =  89.95
xinc         =   0.1
yinc         =  -0.1
xname        = lon
yname        = lat
EOF

# Crear lista ordenada de archivos anuales
ls "$IN_DIR"/SMs_*_GLEAM_v4.2b_MO.nc | sort > "$FILELIST"

echo "Archivos encontrados:"
cat "$FILELIST"

echo "→ Uniendo archivos SMs en un único NetCDF..."

cdo -O -f nc4 -z zip mergetime \
  $(< "$FILELIST") \
  "$MERGED_NC"

echo "→ Remuestreando SMs a 0.25° mediante interpolación bilinear..."

cdo -O -f nc4 -z zip remapbil,"$GRID_025" \
  -setgrid,"$SRC_GRID" \
  "$MERGED_NC" \
  "$OUT_NC"

echo "Proceso finalizado:"
echo "$OUT_NC"



