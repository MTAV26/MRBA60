#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

###############################
# CONFIGURACIÓN
###############################

# Directorio donde están los FWI diarios originales
IN_DIR="/mnt/disco6tb/MHBA60/data/A1_RAW/FWI_v4.1"

# Directorio de salida
OUT_DIR="/mnt/disco6tb/MHBA60/data/A2_TEMP/FWI"

# Fichero de grilla destino 0.25°
GRID_025="/mnt/disco6tb/MHBA60/data/A2_TEMP/grid_0.25.txt"

# Prefijo de los ficheros de entrada
# Formato esperado: fwi-era5-YYYYMMDD.nc
PROD="fwi-era5"

# Periodo de análisis
Y1=2003
Y2=2024

mkdir -p "$OUT_DIR"

cd "$OUT_DIR"

echo "===================================================="
echo "  Flujo FWI ERA5: remapbil → P95 → días > P95"
echo "  Periodo: ${Y1}-${Y2}, producto: ${PROD}"
echo "  Entrada: $IN_DIR"
echo "  Salida:  $OUT_DIR"
echo "===================================================="


###############################
# 0. Comprobaciones
###############################

echo
echo ">> Comprobando entrada y rejilla..."

if [[ ! -f "$GRID_025" ]]; then
  echo "ERROR: no existe el fichero de rejilla:"
  echo "$GRID_025"
  exit 1
fi

test_files=( "$IN_DIR"/${PROD}-2003????.nc )

if [[ ${#test_files[@]} -eq 0 ]]; then
  echo "ERROR: no encuentro ficheros tipo ${PROD}-YYYYMMDD.nc en:"
  echo "$IN_DIR"
  exit 1
fi

printf '%s\n' "${test_files[@]:0:5}"


###############################
# 1. Concatenar diarios por mes en grilla nativa
###############################

echo
echo ">> [1] Concatenando diarios por mes en grilla nativa..."

for m in $(seq -w 1 12); do

  out_native="${PROD}_all_years_month${m}_native.nc"

  if [[ -f "$out_native" ]]; then
    echo "  Mes $m: ya existe ${out_native}, lo dejo."
    continue
  fi

  files_month=( "$IN_DIR"/${PROD}-????${m}??.nc )

  if [[ ${#files_month[@]} -eq 0 ]]; then
    echo "  Mes $m: no hay archivos, salto."
    continue
  fi

  echo "  Mes $m: concatenando ${#files_month[@]} diarios nativos → ${out_native}"

  cdo -s -O -f nc4 -z zip cat \
      "${files_month[@]}" \
      "$out_native"

done


###############################
# 2. Remapeo bilineal a grilla 0.25°
###############################

echo
echo ">> [2] Remapeando a grilla 0.25° con remapbil..."

for m in $(seq -w 1 12); do

  in_native="${PROD}_all_years_month${m}_native.nc"
  out_025="${PROD}_all_years_month${m}_0.25.nc"

  if [[ ! -f "$in_native" ]]; then
    echo "  Mes $m: no existe ${in_native}, salto."
    continue
  fi

  if [[ -f "$out_025" ]]; then
    echo "  Mes $m: ya existe ${out_025}, lo dejo."
    continue
  fi

  echo "  Mes $m: remapbil ${in_native} → ${out_025}"

  cdo -s -O -f nc4 -z zip remapbil,"$GRID_025" \
      "$in_native" \
      "$out_025"

done


###############################
# 3. Percentil 95 mensual en grilla 0.25°
###############################

echo
echo ">> [3] Calculando percentil 95 mensual en 0.25°..."

for m in $(seq -w 1 12); do

  infile="${PROD}_all_years_month${m}_0.25.nc"
  outfile="${PROD}_p95_month${m}_0.25.nc"

  if [[ ! -f "$infile" ]]; then
    echo "  Mes $m: no existe ${infile}, salto."
    continue
  fi

  if [[ -f "$outfile" ]]; then
    echo "  Mes $m: ya existe ${outfile}, lo dejo."
    continue
  fi

  echo "  Mes $m: timpctl,95 sobre ${infile} → ${outfile}"

  cdo -s -O -f nc4 -z zip timpctl,95 \
      "$infile" \
      -timmin "$infile" \
      -timmax "$infile" \
      "$outfile"

done


###############################
# 4. Conteo de días > P95 por año-mes
###############################

echo
echo ">> [4] Contando días con FWI > P95 para cada año-mes..."

mkdir -p tmp_counts

for y in $(seq "$Y1" "$Y2"); do

  echo "  Año $y..."

  for m in $(seq -w 1 12); do

    p95="${PROD}_p95_month${m}_0.25.nc"
    monthfile="${PROD}_all_years_month${m}_0.25.nc"

    tmp_sel="tmp_counts/${PROD}_sel_${y}${m}_0.25.nc"
    tmp_count="tmp_counts/${PROD}_count_${y}${m}_0.25.nc"
    out="tmp_counts/${PROD}_count_${y}${m}_timed_0.25.nc"

    if [[ ! -f "$p95" || ! -f "$monthfile" ]]; then
      echo "    ${y}-${m}: falta ${p95} o ${monthfile}, salto."
      continue
    fi

    if [[ -f "$out" ]]; then
      echo "    ${y}-${m}: ya existe ${out}, lo dejo."
      continue
    fi

    echo "    ${y}-${m}: seleccionando año y contando días > P95..."

    cdo -s -O selyear,"${y}" \
        "$monthfile" \
        "$tmp_sel"

    cdo -s -O timsum -gt \
        "$tmp_sel" \
        "$p95" \
        "$tmp_count"

    cdo -s -O settaxis,${y}-${m}-15,12:00:00,1day \
        "$tmp_count" \
        "$out"

    rm -f "$tmp_sel" "$tmp_count"

  done
done


###############################
# 5. Serie temporal completa
###############################

echo
echo ">> [5] Construyendo serie temporal completa Año-Mes..."

FILELIST="filelist_counts_0.25.txt"
> "$FILELIST"

for y in $(seq "$Y1" "$Y2"); do
  for m in $(seq -w 1 12); do
    f="tmp_counts/${PROD}_count_${y}${m}_timed_0.25.nc"
    if [[ -f "$f" ]]; then
      echo "$f" >> "$FILELIST"
    else
      echo "AVISO: falta $f"
    fi
  done
done

OUT_FINAL="${PROD}_count_exceed95_${Y1}01-${Y2}12_0.25.nc"

echo "  Mergetime sobre $(wc -l < "$FILELIST") ficheros → ${OUT_FINAL}"

cdo -s -O -f nc4 -z zip mergetime \
    $(< "$FILELIST") \
    "$OUT_FINAL"


###############################
# 6. Información rápida
###############################

echo
echo ">> [6] Información rápida del fichero final:"
cdo sinfo "$OUT_FINAL"

echo
echo "===================================================="
echo "  TERMINADO:"
echo "    Fichero final: $OUT_DIR/$OUT_FINAL"
echo "    Contiene nº de días/mes con FWI > P95"
echo "    en rejilla 0.25°, serie ${Y1}-01 ... ${Y2}-12"
echo "===================================================="
