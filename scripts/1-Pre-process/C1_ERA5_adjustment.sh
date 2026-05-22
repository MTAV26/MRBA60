#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Ajuste ERA5 con CDO
# ============================================================

# Rutas
IN_DIR="/mnt/disco6tb/MRBA60/data/A1_RAW"
TMP_DIR="/mnt/disco6tb/MRBA60/data/A2_TEMP"
OUT_DIR="/mnt/disco6tb/MRBA60/data/A3_ADJ"

# Define aquí tu grid objetivo
GRID_DEF="${TMP_DIR}/grid_0.25.txt"

mkdir -p "${TMP_DIR}" "${OUT_DIR}"

# Comprobaciones
command -v cdo >/dev/null 2>&1 || {
  echo "ERROR: cdo no está disponible en el PATH."
  exit 1
}

[[ -f "${GRID_DEF}" ]] || {
  echo "ERROR: no existe el grid file: ${GRID_DEF}"
  exit 1
}

# ============================================================
# Función de procesado
# ============================================================
process_file () {
  local infile="$1"
  local method="$2"

  local inpath="${IN_DIR}/${infile}"

  if [[ ! -f "${inpath}" ]]; then
    echo "AVISO: no existe ${inpath}. Se omite."
    return 0
  fi

  local base="${infile%.nc}"
  local tmpfile="${TMP_DIR}/${base}_TMP.nc"
  local outfile="${OUT_DIR}/${base}_ADJ.nc"

  echo "--------------------------------------------------"
  echo "Procesando: ${infile}"
  echo "Método: ${method}"
  echo "Temporal: ${tmpfile}"
  echo "Salida: ${outfile}"

  # Paso 1: recorte/normalización global
  cdo -L sellonlatbox,-180,180,-90,90 "${inpath}" "${tmpfile}"

  # Paso 2: remapeo al grid objetivo
  case "${method}" in
    remapcon)
      cdo -L remapcon,"${GRID_DEF}" "${tmpfile}" "${outfile}"
      ;;
    remapbil)
      cdo -L remapbil,"${GRID_DEF}" "${tmpfile}" "${outfile}"
      ;;
    *)
      echo "ERROR: método no soportado: ${method}"
      exit 1
      ;;
  esac

  echo "OK: ${outfile}"
}

# ============================================================
# Lista de archivos y método recomendado
# ============================================================
process_file "ERA5-TOT-PREC-2003-2024-MONTLY-025.nc"   "remapcon"
process_file "ERA5-TEMP-MEAN-2003-2024-MONTLY-025.nc"  "remapbil"
process_file "ERA5-WIND-SPEED-2003-2024-MONTLY-025.nc" "remapbil"
process_file "ERA5-TOT-CLOUD-2003-2024-MONTLY-025.nc"  "remapbil"

echo "=================================================="
echo "Proceso completado."
echo "Temporales en: ${TMP_DIR}"
echo "Finales en:    ${OUT_DIR}"
