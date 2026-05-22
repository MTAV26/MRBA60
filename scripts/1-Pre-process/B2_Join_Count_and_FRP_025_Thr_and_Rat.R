rm(list = ls())
graphics.off()
gc()

# ─────────────────────────────────────────────────────────────
# Librerías
# ─────────────────────────────────────────────────────────────
library(terra)

# ─────────────────────────────────────────────────────────────
# Directorios
# ─────────────────────────────────────────────────────────────
in_dir  <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/FRP-MCD14ML-TR"
out_dir <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────
# Periodo objetivo
# ─────────────────────────────────────────────────────────────
start_ym <- "200301"
end_ym   <- "202412"

# ─────────────────────────────────────────────────────────────
# Función para reorganizar arrays
# ─────────────────────────────────────────────────────────────
array_transpuesta <- function(arr) {
  arr <- aperm(arr, c(2, 1, 3))
  arr <- arr[, ncol(arr):1, , drop = FALSE]
  return(arr)
}

# ─────────────────────────────────────────────────────────────
# Función para procesar una variable según conf y ángulo
# ─────────────────────────────────────────────────────────────
procesar_variable <- function(var_prefix, var_name, conf_label, angle_deg,
                              in_dir, out_dir,
                              start_ym = "200301", end_ym = "202412") {
  
  pattern <- paste0(
    "^", var_prefix, "_[0-9]{6}_025_", conf_label,
    "_angle", angle_deg, "deg\\.tif$"
  )
  
  files <- list.files(
    path = in_dir,
    pattern = pattern,
    full.names = TRUE
  )
  
  if (length(files) == 0) {
    message(
      "No se encontraron archivos para ",
      var_name,
      " | ",
      conf_label,
      " | angle ",
      angle_deg
    )
    return(invisible(NULL))
  }
  
  # Extraer YYYYMM
  yyyymm <- sub(
    paste0(
      "^", var_prefix,
      "_([0-9]{6})_025_",
      conf_label,
      "_angle",
      angle_deg,
      "deg\\.tif$"
    ),
    "\\1",
    basename(files)
  )
  
  # Filtrar directamente al periodo seleccionado
  keep <- yyyymm >= start_ym & yyyymm <= end_ym
  files <- files[keep]
  yyyymm <- yyyymm[keep]
  
  if (length(files) == 0) {
    message(
      "⚠️ No hay archivos dentro del periodo ",
      start_ym,
      "-",
      end_ym,
      " para ",
      var_name,
      " | ",
      conf_label,
      " | angle ",
      angle_deg
    )
    return(invisible(NULL))
  }
  
  # Ordenar
  ord <- order(yyyymm)
  files <- files[ord]
  yyyymm <- yyyymm[ord]
  
  message(
    "Procesando: ",
    var_name,
    " | ",
    conf_label,
    " | angle ",
    angle_deg
  )
  
  cat("Fechas detectadas:\n", paste(yyyymm, collapse = ", "), "\n")
  
  # Serie mensual esperada
  dates_full <- seq(
    as.Date(paste0(start_ym, "01"), format = "%Y%m%d"),
    as.Date(paste0(end_ym,   "01"), format = "%Y%m%d"),
    by = "month"
  )
  
  yyyymm_full <- format(dates_full, "%Y%m")
  
  # Comprobar meses ausentes
  missing_ym <- setdiff(yyyymm_full, yyyymm)
  
  if (length(missing_ym) > 0) {
    warning(
      "Faltan meses en ",
      var_name,
      " | ",
      conf_label,
      " | angle ",
      angle_deg,
      ": ",
      paste(missing_ym, collapse = ", ")
    )
  }
  
  # Leer stack y convertir a array
  r_stack <- rast(files)
  r_array <- as.array(r_stack)
  r_array <- array_transpuesta(r_array)
  
  # Nombre del objeto
  output_name <- paste0(
    var_name,
    "_",
    conf_label,
    "_angle",
    angle_deg
  )
  
  # Guardar en entorno global
  assign(output_name, r_array, envir = .GlobalEnv)
  
  out_file <- file.path(
    out_dir,
    paste0(
      "MODIS-",
      output_name,
      "-",
      start_ym,
      "-",
      end_ym,
      "-025.RData"
    )
  )
  
  save(list = output_name, file = out_file)
  
  message("Guardado: ", out_file)
  
  rm(
    r_stack,
    r_array,
    files,
    yyyymm,
    dates_full,
    yyyymm_full,
    missing_ym,
    keep,
    ord
  )
  gc()
  
  return(invisible(NULL))
}

# ─────────────────────────────────────────────────────────────
# Parámetros
# ─────────────────────────────────────────────────────────────
conf_labels  <- "conf30"
angle_levels <- 30

# Si más adelante quieres varios:
# conf_levels <- c(0, seq(50, 100, by = 5))
# conf_labels <- ifelse(conf_levels == 0, "nocut", paste0("conf", conf_levels))
# angle_levels <- 30

# ─────────────────────────────────────────────────────────────
# Procesamiento automático
# ─────────────────────────────────────────────────────────────
for (conf_label in conf_labels) {
  for (angle_deg in angle_levels) {
    
    # ---------------------------------------------------------
    # FRP suma mensual
    # ---------------------------------------------------------
    procesar_variable(
      var_prefix = "MCD14ML_FRPsum",
      var_name   = "FRPsum",
      conf_label = conf_label,
      angle_deg  = angle_deg,
      in_dir     = in_dir,
      out_dir    = out_dir,
      start_ym   = start_ym,
      end_ym     = end_ym
    )
    
    # ---------------------------------------------------------
    # FRP media mensual
    # ---------------------------------------------------------
    procesar_variable(
      var_prefix = "MCD14ML_FRPmean",
      var_name   = "FRPmean",
      conf_label = conf_label,
      angle_deg  = angle_deg,
      in_dir     = in_dir,
      out_dir    = out_dir,
      start_ym   = start_ym,
      end_ym     = end_ym
    )
    
    # ---------------------------------------------------------
    # FRP mediana mensual
    # ---------------------------------------------------------
    procesar_variable(
      var_prefix = "MCD14ML_FRPmedian",
      var_name   = "FRPmedian",
      conf_label = conf_label,
      angle_deg  = angle_deg,
      in_dir     = in_dir,
      out_dir    = out_dir,
      start_ym   = start_ym,
      end_ym     = end_ym
    )
    
    # ---------------------------------------------------------
    # Recuento mensual de Active Fires
    # ---------------------------------------------------------
    procesar_variable(
      var_prefix = "MCD14ML_AFcount",
      var_name   = "AFcount",
      conf_label = conf_label,
      angle_deg  = angle_deg,
      in_dir     = in_dir,
      out_dir    = out_dir,
      start_ym   = start_ym,
      end_ym     = end_ym
    )
  }
}
