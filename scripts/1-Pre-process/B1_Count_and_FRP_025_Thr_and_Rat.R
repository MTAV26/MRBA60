rm(list = ls())
gc()

# =========================================================
# LIBRERÍAS
# =========================================================
library(data.table)
library(terra)

# =========================================================
# RUTAS
# =========================================================
in_dir  <- "/mnt/disco6tb/MRBA60/data/A1_RAW/MCD14ML"
out_dir <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/FRP-MCD14ML-TR"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# GRID GLOBAL 0.25°
# =========================================================
grid_raster_global_template <- rast(
  ncols = 1440, nrows = 720,
  xmin = -180, xmax = 180,
  ymin = -90, ymax = 90,
  crs = "EPSG:4326"
)

# =========================================================
# PARÁMETROS
# =========================================================
# Constante para convertir sample a ángulo de escaneo
s <- 0.0014184397

# Umbrales de confianza
# conf_levels <- c(0, seq(50, 100, by = 5))
conf_levels <- 30

# Umbral(es) de ángulo en grados
angle_levels_deg <- 30
angle_levels_rad <- angle_levels_deg * pi / 180

# =========================================================
# LISTAR ARCHIVOS DISPONIBLES
# =========================================================
# Acepta tanto 006 como 061
files_in <- list.files(
  path = in_dir,
  pattern = "^MCD14ML\\.[0-9]{6}\\.(006|061)\\.03\\.txt\\.gz$",
  full.names = TRUE
)

if (length(files_in) == 0) {
  stop("No se encontraron archivos MCD14ML en: ", in_dir)
}

files_in <- sort(files_in)

# Extraer YYYYMM del nombre
yyyymm <- sub(
  "^MCD14ML\\.([0-9]{6})\\.(006|061)\\.03\\.txt\\.gz$",
  "\\1",
  basename(files_in)
)

# =========================================================
# PROCESAMIENTO MENSUAL
# =========================================================
for (k in seq_along(files_in)) {
  
  file_path <- files_in[k]
  ym <- yyyymm[k]
  
  message("Procesando: ", basename(file_path))
  
  # ---------------------------------
  # Leer datos
  # ---------------------------------
  data <- fread(
    cmd = paste("zcat", shQuote(file_path)),
    skip = 1,
    header = FALSE
  )
  
  colnames(data) <- c(
    "YYYYMMDD", "HHMM", "sat", "lat", "lon", "T21", "T31",
    "sample", "FRP", "conf", "type", "dn"
  )
  
  num_cols <- c("lat", "lon", "T21", "T31", "sample", "FRP", "conf")
  data[, (num_cols) := lapply(.SD, as.numeric), .SDcols = num_cols]
  
  # ---------------------------------
  # Filtrar solo fuegos tipo 0
  # ---------------------------------
  data <- data[type == 0]
  
  if (nrow(data) == 0) {
    message("Sin datos tipo == 0 en ", ym)
    rm(data)
    gc()
    next
  }
  
  # ---------------------------------
  # Calcular ángulo de escaneo
  # ---------------------------------
  data[, angle_rad := (sample - 676.5) * s]
  data[, angle_abs_rad := abs(angle_rad)]
  
  # =====================================================
  # BUCLE POR CONFIANZA
  # =====================================================
  for (conf_threshold in conf_levels) {
    
    data_conf <- data[conf >= conf_threshold]
    
    if (nrow(data_conf) == 0) {
      message("Sin datos para conf >= ", conf_threshold, " en ", ym)
      next
    }
    
    # =================================================
    # BUCLE POR ÁNGULO
    # =================================================
    for (i in seq_along(angle_levels_rad)) {
      
      angle_threshold <- angle_levels_rad[i]
      ang_deg <- angle_levels_deg[i]
      
      data_sub <- data_conf[angle_abs_rad <= angle_threshold]
      
      if (nrow(data_sub) == 0) {
        message(
          "Sin datos para angle <= ", ang_deg,
          "° y conf >= ", conf_threshold,
          " en ", ym
        )
        next
      }
      
      # ---------------------------------
      # Añadir campo contador de AF
      # ---------------------------------
      data_sub[, AF_count := 1L]
      
      # ---------------------------------
      # Convertir a vector espacial
      # ---------------------------------
      data_vect <- vect(
        data_sub,
        geom = c("lon", "lat"),
        crs = "EPSG:4326"
      )
      
      # ---------------------------------
      # Raster FRP sum
      # ---------------------------------
      grid_raster_sum <- rasterize(
        x = data_vect,
        y = grid_raster_global_template,
        field = "FRP",
        fun = "sum",
        background = NA
      )
      
      # ---------------------------------
      # Raster FRP mean
      # ---------------------------------
      grid_raster_mean <- rasterize(
        x = data_vect,
        y = grid_raster_global_template,
        field = "FRP",
        fun = "mean",
        background = NA
      )
      
      # ---------------------------------
      # Raster FRP median
      # ---------------------------------
      grid_raster_median <- rasterize(
        x = data_vect,
        y = grid_raster_global_template,
        field = "FRP",
        fun = "median",
        background = NA
      )
      
      # ---------------------------------
      # Raster recuento de Active Fires
      # ---------------------------------
      grid_raster_count <- rasterize(
        x = data_vect,
        y = grid_raster_global_template,
        field = "AF_count",
        fun = "sum",
        background = NA
      )
      
      # ---------------------------------
      # Etiqueta de confianza
      # ---------------------------------
      conf_label <- ifelse(
        conf_threshold == 0,
        "nocut",
        paste0("conf", conf_threshold)
      )
      
      # ---------------------------------
      # Nombres de salida
      # ---------------------------------
      out_file_sum <- file.path(
        out_dir,
        paste0(
          "MCD14ML_FRPsum_",
          ym,
          "_025_",
          conf_label,
          "_angle",
          ang_deg,
          "deg.tif"
        )
      )
      
      out_file_mean <- file.path(
        out_dir,
        paste0(
          "MCD14ML_FRPmean_",
          ym,
          "_025_",
          conf_label,
          "_angle",
          ang_deg,
          "deg.tif"
        )
      )
      
      out_file_median <- file.path(
        out_dir,
        paste0(
          "MCD14ML_FRPmedian_",
          ym,
          "_025_",
          conf_label,
          "_angle",
          ang_deg,
          "deg.tif"
        )
      )
      
      out_file_count <- file.path(
        out_dir,
        paste0(
          "MCD14ML_AFcount_",
          ym,
          "_025_",
          conf_label,
          "_angle",
          ang_deg,
          "deg.tif"
        )
      )
      
      # ---------------------------------
      # Guardar
      # ---------------------------------
      writeRaster(
        grid_raster_sum,
        filename = out_file_sum,
        overwrite = TRUE
      )
      
      writeRaster(
        grid_raster_mean,
        filename = out_file_mean,
        overwrite = TRUE
      )
      
      writeRaster(
        grid_raster_median,
        filename = out_file_median,
        overwrite = TRUE
      )
      
      writeRaster(
        grid_raster_count,
        filename = out_file_count,
        overwrite = TRUE
      )
      
      message(
        "Guardados: ",
        ym,
        " | ",
        conf_label,
        " | angle <= ",
        ang_deg,
        "°"
      )
      
      # ---------------------------------
      # Limpiar
      # ---------------------------------
      rm(
        data_vect,
        grid_raster_sum,
        grid_raster_mean,
        grid_raster_median,
        grid_raster_count
      )
      gc()
    }
    
    rm(data_conf)
    gc()
  }
  
  rm(data)
  gc()
}
