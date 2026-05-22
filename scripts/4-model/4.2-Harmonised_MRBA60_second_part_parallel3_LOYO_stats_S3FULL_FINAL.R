# ==============================================================================
# MRBA60 / FireCCI60
# SEGUNDA PARTE DESDE CSV DE CHECK:
#   1) Lee SelectedPredictors_<bioma>_COMMON.csv
#   2) Mantiene LOYO comentado en el código
#   3) Entrena modelo común final por mes
#   4) Estima SHAP por bioma y mes
#   5) Aplica modelo completo retrospectivo 2003-2024
#   6) Procesa biomas en paralelo y une los mosaicos globales al final
#
# NO recalcula autocorrelación ni RFE.
# GPP no se carga ni entra nunca en el modelo.
# ==============================================================================

rm(list = ls())
graphics.off()
gc()

# Evitar sobre-paralelización interna de BLAS/OpenMP
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

# ==============================================================================
# 0. PAQUETES
# ==============================================================================

required_packages <- c(
  "dplyr",
  "lubridate",
  "sf",
  "ncdf4",
  "randomForest",
  "parallel",
  "fastshap",
  "tidyr",
  "tibble",
  "ggplot2",
  "RColorBrewer",
  "scales",
  "grid"
)

missing_packages <- required_packages[
  !sapply(required_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Faltan paquetes: ",
    paste(missing_packages, collapse = ", "),
    "\nInstálalos antes de ejecutar el script."
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(sf)
  library(ncdf4)
  library(randomForest)
  library(parallel)
  library(fastshap)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(grid)
})

sf::sf_use_s2(FALSE)

# ==============================================================================
# 1. CONFIGURACIÓN GENERAL
# ==============================================================================

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

output_dir          <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv      <- file.path(output_dir, "csv")
output_dir_plot     <- file.path(output_dir, "plot")
output_dir_plot_rle <- file.path(output_dir, "plot_rle")
output_dir_plot_sc  <- file.path(output_dir, "plot_scatter")
output_dir_RData    <- file.path(output_dir, "RData")
output_dir_log      <- file.path(output_dir, "logs_harmonisation")

selected_predictors_dir <- output_dir_csv
selected_predictors_dir2 = "/mnt/disco6tb/MRBA60/results/B1-MRBA60-2003-2024/csv-1/"
dirs <- c(
  output_dir,
  output_dir_csv,
  output_dir_plot,
  output_dir_plot_rle,
  output_dir_plot_sc,
  output_dir_RData,
  output_dir_log
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Número de núcleos: un bioma por núcleo.
# Los mosaicos globales se reconstruyen después, en una fase secuencial.
N_CORES <- min(4, max(1, parallel::detectCores() - 1))

SEED <- 123
N_SHAP <- 100
SAVE_TIMESERIES_PLOTS <- TRUE
SAVE_MODELS <- FALSE

cat("N_CORES =", N_CORES, "\n")

# Biomas a procesar.
# Usa NULL para todos, o vector numérico para solo algunos.
# idx_biomas_pendientes <- unique(c(7, 8)) #NULL
# idx_biomas_pendientes <- unique(c(17:50, 9, 7, 8)) #NULL
# idx_biomas_pendientes <- c(setdiff(1:50, c(7, 8)), 7, 8)
idx_biomas_pendientes <- NULL
# idx_biomas_pendientes <- c(22,23)
# ==============================================================================
# 2. FUNCIONES AUXILIARES
# ==============================================================================
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_nc_var.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/safe_name.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/apply_mask_3d_2.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/build_prediction_df.R")

source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/parse_predictor_string.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/load_selected_predictors.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/calc_stats.R")

source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/format_stat_value.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/format_stats_text.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/format_full_summary_text.R")

source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/predict_layer_rf.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/build_monthly_training_data.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/build_monthly_training_data.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/plot_common_timeseries.R")

# ==============================================================================
# 3. CARGA DE DATOS GLOBALES
# ==============================================================================

cat("Cargando lon/lat...\n")
load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))

dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

cat("Cargando biomas...\n")
biomas_shp <- sf::st_read(
  file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"),
  quiet = TRUE
)

biomas_shp <- sf::st_transform(biomas_shp, crs = 4326)

biomas_all <- unique(biomas_shp$cont_bm)

if (is.null(idx_biomas_pendientes)) {
  biomas_unique <- biomas_all
} else {
  biomas_unique <- biomas_all[idx_biomas_pendientes]
}

cat("Biomas a procesar:", length(biomas_unique), "\n")

# ------------------------------------------------------------------------------
# Active Fire count
# ------------------------------------------------------------------------------

cat("Cargando Active Fire count...\n")
load(file.path(dir_oss, "MODIS-AFcount_conf30_angle30-200301-202412-025.RData"))

count_ActiveFire_tot <- AFcount_conf30_angle30
count_ActiveFire_tot[count_ActiveFire_tot == 0] <- NA

rm(AFcount_conf30_angle30)
gc()

# ------------------------------------------------------------------------------
# FireCCIS311
# ------------------------------------------------------------------------------

cat("Cargando FireCCIS311...\n")
load(file.path(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))

BA_FireS3 <- s3 / 1e6
BA_FireS3[BA_FireS3 == 0] <- NA

rm(s3)
gc()

# ------------------------------------------------------------------------------
# FireCCI51
# ------------------------------------------------------------------------------

cat("Cargando FireCCI51...\n")
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

BA_Fire51_tot <- f51 / 1e6
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA

rm(f51)
gc()

# ------------------------------------------------------------------------------
# FWI95d
# ------------------------------------------------------------------------------

cat("Cargando FWI95d...\n")
load(file.path(dir_oss, "FWI95d-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"))

FWI_tot <- FWI95d_tot
rm(FWI95d_tot)
gc()

# ------------------------------------------------------------------------------
# ERA5 TEMP
# ------------------------------------------------------------------------------

cat("Cargando temperatura...\n")
TEMP_tot <- read_nc_var(
  file.path(dir_oss, "ERA5-TEMP-MEAN-2003-2024-MONTLY-025_ADJ.nc"),
  "t2m"
)

TEMP_tot <- TEMP_tot - 273.15
gc()

# ------------------------------------------------------------------------------
# ERA5 PREC
# ------------------------------------------------------------------------------

cat("Cargando precipitación...\n")
PREC_tot <- read_nc_var(
  file.path(dir_oss, "ERA5-TOT-PREC-2003-2024-MONTLY-025_ADJ.nc"),
  "tp"
)

PREC_tot <- PREC_tot * 1000
gc()

# ------------------------------------------------------------------------------
# ERA5 WIND
# ------------------------------------------------------------------------------

cat("Cargando viento...\n")
WIND_tot <- read_nc_var(
  file.path(dir_oss, "ERA5-WIND-SPEED-2003-2024-MONTLY-025_ADJ.nc"),
  "si10"
)

gc()

# ------------------------------------------------------------------------------
# NDVI
# ------------------------------------------------------------------------------

cat("Cargando NDVI...\n")
load(file.path(dir_oss, "NDVI-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"))
gc()

# ------------------------------------------------------------------------------
# FRP median / sum
# ------------------------------------------------------------------------------

cat("Cargando FRP median...\n")
load(file.path(dir_oss, "MODIS-FRPmedian_conf30_angle30-200301-202412-025.RData"))

FRPmedian_tot <- FRPmedian_conf30_angle30
rm(FRPmedian_conf30_angle30)
gc()

cat("Cargando FRP sum...\n")
load(file.path(dir_oss, "MODIS-FRPsum_conf30_angle30-200301-202412-025.RData"))

FRPsum_tot <- FRPsum_conf30_angle30
rm(FRPsum_conf30_angle30)
gc()

# ------------------------------------------------------------------------------
# ERA5 cloud
# ------------------------------------------------------------------------------

cat("Cargando cloud cover...\n")
cloud_tot <- read_nc_var(
  file.path(dir_oss, "ERA5-TOT-CLOUD-2003-2024-MONTLY-025_ADJ.nc"),
  "tcc"
)

gc()

# ------------------------------------------------------------------------------
# VPD
# ------------------------------------------------------------------------------

cat("Cargando VPD...\n")
load(file.path(dir_oss, "VPD-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"))
gc()

# ------------------------------------------------------------------------------
# Soil moisture
# ------------------------------------------------------------------------------

cat("Cargando soil moisture...\n")
load(file.path(dir_oss, "SMs-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"))

SOIL_tot <- SMs_tot
rm(SMs_tot)
gc()

# ==============================================================================
# 4. RECONSTRUIR FIRECCIS311 EN ARRAY COMPLETO 2003-2024
# ==============================================================================

nrows <- dim(BA_Fire51_tot)[1]
ncols <- dim(BA_Fire51_tot)[2]
ntime_full <- length(dates_full)

if (length(ind_common) != dim(BA_FireS3)[3]) {
  stop(
    "La longitud de ind_common no coincide con la dimensión temporal de BA_FireS3."
  )
}

BA_FireS3_tot <- array(
  NA_real_,
  dim = c(nrows, ncols, ntime_full)
)

BA_FireS3_tot[, , ind_common] <- BA_FireS3

rm(BA_FireS3)
gc()

# ==============================================================================
# 5. COORDENADAS Y ÁREA DE CELDA
# ==============================================================================

nlon <- length(lon)
nlat <- length(lat)

lon_mat <- matrix(
  rep(lon, each = nlat),
  nrow = nlat,
  ncol = nlon,
  byrow = FALSE
)

lat_mat <- matrix(
  rep(lat, times = nlon),
  nrow = nlat,
  ncol = nlon,
  byrow = FALSE
)

lon_range <- lon_mat[1, ]
lat_range <- lat_mat[, 1]

cell_area_constant <- (110.57 * 0.25) * (111.32 * 0.25)
area_by_row <- cell_area_constant * cos(lat * pi / 180)

area_matrix <- matrix(
  rep(area_by_row, times = nlon),
  nrow = nlat,
  ncol = nlon,
  byrow = FALSE
)

area_matrix <- t(area_matrix)

# ==============================================================================
# 6. CARGAR MÁSCARA DE FUEGO Y AJUSTAR CEROS INTERNOS
# ==============================================================================

cat("Cargando FireMask_AF3030F...\n")

dir_mask <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData"
file_firemask <- file.path(dir_mask, "FireMask_AF3030F.RData")

load(file_firemask)

status_matrix2_tot <- FireMask_AF3030F
rm(FireMask_AF3030F)
gc()

inside <- (!is.na(status_matrix2_tot)) & (status_matrix2_tot == 1)

BA_Fire51_tot[inside & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[inside & is.na(BA_FireS3_tot)] <- 0

count_ActiveFire_tot[!inside] <- NA
count_ActiveFire_tot[inside & is.na(count_ActiveFire_tot)] <- 0

FRPmedian_tot[!inside] <- NA
FRPmedian_tot[inside & is.na(FRPmedian_tot)] <- 0

FRPsum_tot[!inside] <- NA
FRPsum_tot[inside & is.na(FRPsum_tot)] <- 0

rm(inside)
gc()

# ==============================================================================
# 7. PREDICTORES COMUNES Y ARRAYS GLOBALES
# ==============================================================================

PREDICTOR_LABELS <- c(
  "f5" = "FireCCI51",
  "count_ActiveFire" = "Active Fire",
  "prec" = "Precipitation",
  "temp" = "Temperature",
  "FRPsum" = "FRP sum",
  "FRPmedian" = "FRP median",
  "NDVI" = "NDVI",
  "FWI" = "FWI",
  "wind" = "Wind Speed",
  "lat" = "Latitude",
  "lon" = "Longitude",
  "cloud" = "Cloud",
  "vpd" = "VPD",
  "soil" = "Soil moisture"
)

candidate_predictors <- names(PREDICTOR_LABELS)
response <- "f3"

# Los arrays globales se inicializan después de la fase paralela.
# Cada proceso hijo guarda sus arrays locales por bioma y el proceso padre
# los une secuencialmente para evitar escrituras concurrentes y copias innecesarias.
#
global_BA_FireHarmonized_common_loyo <- array(
  NA_real_,
  dim = c(nrows, ncols, length(ind_common))
)

# ==============================================================================
# 8. FUNCIÓN PRINCIPAL: ARMONIZACIÓN DESDE CSV DE PREDICTORES
# ==============================================================================

process_biome_harmonisation <- function(bioma) {
  
  t_start <- Sys.time()
  safe_biome_name <- safe_name(as.character(bioma))
  
  log_file <- file.path(output_dir_log, paste0("LOG_Harmonisation_", safe_biome_name, ".txt"))
  log_con <- file(log_file, open = "wt")
  
  sink(log_con)
  sink(log_con, type = "message")
  
  on.exit({
    sink(type = "message")
    sink()
    close(log_con)
  }, add = TRUE)
  
  cat("============================================================\n")
  cat("Procesando armonización para bioma:", bioma, "\n")
  cat("Inicio:", as.character(t_start), "\n")
  cat("PID:", Sys.getpid(), "\n")
  cat("============================================================\n")
  
  selected_info <- load_selected_predictors(safe_biome_name)
  
  if (is.null(selected_info)) {
    cat("No existe CSV de predictores para el bioma:", safe_biome_name, "\n")
    
    return(list(
      Biome = safe_biome_name,
      Status = "SKIPPED_NO_SELECTED_PREDICTORS_CSV",
      n_months_modelled = 0L,
      seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    ))
  }
  
  cat("CSV de predictores leído:\n")
  cat(selected_info$file, "\n")
  
  # ---------------------------------------------------------------------------
  # 8.1 Selección espacial del bioma
  # ---------------------------------------------------------------------------
  
  bioma_sel <- biomas_shp %>%
    dplyr::filter(cont_bm == bioma)
  
  bbox <- sf::st_bbox(bioma_sel)
  
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  
  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    cat("Sin celdas para el bioma:", bioma, "\n")
    
    return(list(
      Biome = safe_biome_name,
      Status = "SKIPPED_NO_GRID_CELLS",
      n_months_modelled = 0L,
      seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    ))
  }
  
  crop3 <- function(x) {
    x[lon_idx, lat_idx, , drop = FALSE]
  }
  
  status_crop <- crop3(status_matrix2_tot)
  
  lon_mat_crop <- lon_mat[lat_idx, lon_idx, drop = FALSE]
  lat_mat_crop <- lat_mat[lat_idx, lon_idx, drop = FALSE]
  
  lon_vec <- lon[lon_idx]
  lat_vec <- lat[lat_idx]
  lon_vec_crop <- lon_vec
  lat_vec_crop <- lat_vec
  
  # ---------------------------------------------------------------------------
  # 8.2 Asignación del bioma a la grilla
  # ---------------------------------------------------------------------------
  
  grid_points_crop <- sf::st_as_sf(
    data.frame(
      lon = as.vector(lon_mat_crop),
      lat = as.vector(lat_mat_crop)
    ),
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  inter <- sf::st_intersects(grid_points_crop, biomas_shp)
  
  bioma_asignado <- sapply(inter, function(i) {
    if (length(i) == 0) return("Ninguno")
    biomas_shp$cont_bm[i[1]]
  })
  
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  presencia_en_serie <- apply(
    status_crop,
    c(1, 2),
    function(x) any(x == 1, na.rm = TRUE)
  )
  
  grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie))
  
  puntos_sin_bioma_con_presencia <- grid_points_biomas %>%
    dplyr::filter(bioma_final == "Ninguno" & presencia_serie)
  
  puntos_con_bioma <- grid_points_biomas %>%
    dplyr::filter(bioma_final != "Ninguno")
  
  if (
    nrow(puntos_sin_bioma_con_presencia) > 0 &&
      nrow(puntos_con_bioma) > 0
  ) {
    
    nearest_idx <- sf::st_nearest_feature(
      puntos_sin_bioma_con_presencia,
      puntos_con_bioma
    )
    
    grid_points_biomas$bioma_final[
      which(
        grid_points_biomas$bioma_final == "Ninguno" &
          grid_points_biomas$presencia_serie
      )
    ] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  mask <- grid_points_biomas$bioma_final == bioma
  
  mask_matrix <- matrix(
    mask,
    nrow = length(lat_vec),
    ncol = length(lon_vec),
    byrow = FALSE
  )
  
  mask_lonlat <- t(mask_matrix)
  
  if (!any(mask_lonlat, na.rm = TRUE)) {
    cat("Máscara vacía para bioma:", bioma, "\n")
    
    return(list(
      Biome = safe_biome_name,
      Status = "SKIPPED_EMPTY_MASK",
      n_months_modelled = 0L,
      seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    ))
  }
  
  # ---------------------------------------------------------------------------
  # 8.3 Recorte y máscara de predictores
  # ---------------------------------------------------------------------------
  
  cat("Recortando arrays...\n")
  
  arr <- list(
    s3 = crop3(BA_FireS3_tot),
    f51 = crop3(BA_Fire51_tot),
    af = crop3(count_ActiveFire_tot),
    prec = crop3(PREC_tot),
    temp = crop3(TEMP_tot),
    frp_sum = crop3(FRPsum_tot),
    frp_median = crop3(FRPmedian_tot),
    fwi = crop3(FWI_tot),
    ndvi = crop3(NDVI_tot),
    wind = crop3(WIND_tot),
    cloud = crop3(cloud_tot),
    vpd = crop3(VPD_tot),
    soil = crop3(SOIL_tot)
  )
  
  arr <- lapply(arr, apply_mask_3d_2, mask_lonlat = mask_lonlat)
  
  lat_2d <- t(lat_mat_crop)
  lon_2d <- t(lon_mat_crop)
  
  lat_2d[!mask_lonlat] <- NA
  lon_2d[!mask_lonlat] <- NA
  
  arr$lat <- array(
    rep(lat_2d, ntime_full),
    dim = c(dim(lat_2d), ntime_full)
  )
  
  arr$lon <- array(
    rep(lon_2d, ntime_full),
    dim = c(dim(lon_2d), ntime_full)
  )
  
  arr_common <- lapply(
    arr,
    function(x) x[, , ind_common, drop = FALSE]
  )
  
  area_matrix_crop <- area_matrix[lon_idx, lat_idx, drop = FALSE]
  area_matrix_crop[!mask_lonlat] <- NA
  
  total_fire <- sum(arr_common$s3, na.rm = TRUE)
  
  if (!is.finite(total_fire) || total_fire == 0) {
    cat("El bioma no registra incendios en el periodo common. Se omite.\n")
    
    return(list(
      Biome = safe_biome_name,
      Status = "SKIPPED_NO_COMMON_FIRE",
      n_months_modelled = 0L,
      seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    ))
  }
  
  # ---------------------------------------------------------------------------
  # 8.4 Inicialización de modelos, tablas y salidas locales
  # ---------------------------------------------------------------------------
  
  time_seq_common <- dates_common
  time_seq_full <- dates_full
  ntime_common <- length(time_seq_common)
  
  obs_BA_FireS3_biome <- arr_common$s3
  obs_BA_Fire51_biome <- arr_common$f51
  
  obs_BA_FireS3_biome_full <- arr$s3
  obs_BA_Fire51_biome_full <- arr$f51
  
   BA_FireHarmonized_common_loyo <- obs_BA_Fire51_biome
  BA_FireHarmonized_common <- obs_BA_Fire51_biome
  BA_FireHarmonized_full <- obs_BA_Fire51_biome_full
  
  final_models <- vector("list", 12)
  final_predictors <- vector("list", 12)
  
  evaluation_log_loyo <- data.frame(
    Month = integer(),
    Model = character(),
    Formula = character(),
    RMSE = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  evaluation_log <- data.frame(
    Month = integer(),
    Model = character(),
    Formula = character(),
    RMSE = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )
  
  table_selected_used <- data.frame(
    Bioma = character(),
    Month = integer(),
    Predictor_Source_CSV = character(),
    Predictors_used = character(),
    Formula_used = character(),
    n_predictors_used = integer(),
    n_training_rows = integer(),
    n_valid_years_ge30 = integer(),
    stringsAsFactors = FALSE
  )
  
  df_shap_global <- list()
  
  f5_max_train_by_month <- rep(NA_real_, 12)
  f3_max_train_by_month <- rep(NA_real_, 12)
  tope_max_train_by_month <- rep(NA_real_, 12)
  
  for (m in 1:12) {
    idx_m <- which(lubridate::month(time_seq_common) == m)
    
    if (length(idx_m) > 0) {
      v_f5 <- as.vector(obs_BA_Fire51_biome[, , idx_m])
      v_f3 <- as.vector(obs_BA_FireS3_biome[, , idx_m])
      
      f5_max_train_by_month[m] <- suppressWarnings(max(v_f5, na.rm = TRUE))
      f3_max_train_by_month[m] <- suppressWarnings(max(v_f3, na.rm = TRUE))
      tope_max_train_by_month[m] <- suppressWarnings(max(c(v_f5, v_f3), na.rm = TRUE))
    }
  }
  
  max_train_tbl <- data.frame(
    Biome = safe_biome_name,
    Month = 1:12,
    Max_F51_in_COMMON = f5_max_train_by_month,
    Max_F3_in_COMMON = f3_max_train_by_month,
    Tope_MAX_F51_F3 = tope_max_train_by_month
  )
  
  write.csv(
    max_train_tbl,
    file = file.path(output_dir_csv, paste0("Max_Tope_COMMON_ByMonth_", safe_biome_name, ".csv")),
    row.names = FALSE
  )
  
  # ---------------------------------------------------------------------------
  # 8.5 Modelo común y SHAP por mes leyendo predictores del CSV. LOYO queda comentado.
  # ---------------------------------------------------------------------------
  
  n_months_modelled <- 0L
  
  for (mes in 1:12) {
    
    cat("\nProcesando mes desde CSV de predictores:", mes, "\n")
    
    best_preds <- selected_info$by_month[[mes]]
    best_preds <- unique(best_preds[!tolower(best_preds) %in% "gpp"])
    best_preds <- intersect(best_preds, candidate_predictors)
    
    month_data <- build_monthly_training_data(
      arr_common = arr_common,
      mes = mes
    )
    
    mes_indices_central <- month_data$mes_indices_central
    df_full <- month_data$df_full
    
    if (length(best_preds) == 0) {
      cat("Mes", mes, ": no hay predictores RFE en el CSV. Se conserva F51.\n")
      next
    }
    
    if (is.null(df_full) || nrow(df_full) == 0) {
      cat("Mes", mes, ": ventana sin datos. Se conserva F51.\n")
      next
    }
    
    best_preds <- best_preds[
      best_preds %in% names(df_full)
    ]
    
    if (length(best_preds) == 0) {
      cat("Mes", mes, ": predictores del CSV no están en df_full. Se conserva F51.\n")
      next
    }
    
    df_full <- df_full[
      complete.cases(df_full[, c(response, best_preds), drop = FALSE]),
    ]
    
    if (nrow(df_full) == 0) {
      cat("Mes", mes, ": sin casos completos para los predictores del CSV. Se conserva F51.\n")
      next
    }
    
    var_check <- sapply(
      df_full[, best_preds, drop = FALSE],
      function(x) stats::var(x, na.rm = TRUE)
    )
    
    best_preds <- best_preds[
      is.finite(var_check) & var_check > 0
    ]
    
    if (length(best_preds) == 0) {
      cat("Mes", mes, ": ningún predictor del CSV con varianza > 0. Se conserva F51.\n")
      next
    }
    
    year_counts <- table(df_full$year)
    n_valid_years_ge30 <- sum(year_counts >= 30)
    
    if (n_valid_years_ge30 < 2) {
      cat("Mes", mes, ": datos insuficientes por año. Se conserva F51.\n")
      next
    }
    
    formula_modelo <- as.formula(
      paste(response, "~", paste(best_preds, collapse = " + "))
    )
    
    mtry_value <- min(length(best_preds), max(1, round(length(best_preds) / 3)))
    
    table_selected_used <- rbind(
      table_selected_used,
      data.frame(
        Bioma = safe_biome_name,
        Month = mes,
        Predictor_Source_CSV = selected_info$file,
        Predictors_used = paste(best_preds, collapse = " + "),
        Formula_used = deparse(formula_modelo),
        n_predictors_used = length(best_preds),
        n_training_rows = nrow(df_full),
        n_valid_years_ge30 = n_valid_years_ge30,
        stringsAsFactors = FALSE
      )
    )
    
    cat("Mes", mes, "\n")
    cat("  Predictores usados desde CSV:", paste(best_preds, collapse = " + "), "\n")
    cat("  Filas entrenamiento:", nrow(df_full), "\n")
    cat("  Años válidos >=30:", n_valid_years_ge30, "\n")
    
#     # ---------------------- LOYO ----------------------
    for (target_year in unique(df_full$year)) {

      cat("  - Excluyendo año:", target_year, "para LOYO\n")

      train_data <- subset(df_full, year != target_year)
      train_data <- train_data[
        complete.cases(train_data[, c(response, best_preds), drop = FALSE]),
      ]

      if (nrow(train_data) < 30 || length(unique(train_data$year)) < 2) {
        cat("    Insuficiente tras excluir", target_year, ".\n")
        next
      }

      set.seed(SEED + mes * 1000 + as.integer(target_year))
      modelo_loyo <- randomForest::randomForest(
        formula_modelo,
        data = train_data,
        ntree = 250,
        mtry = mtry_value,
        importance = TRUE
      )

      idx_time <- which(lubridate::year(time_seq_common[mes_indices_central]) == target_year)

      for (ii in idx_time) {
        tt <- mes_indices_central[ii]

        BA_FireHarmonized_common_loyo[, , tt] <- predict_layer_rf(
          modelo = modelo_loyo,
          arr = arr_common,
          tt = tt,
          best_preds = best_preds,
          area_matrix_crop = area_matrix_crop,
          base_layer = BA_FireHarmonized_common_loyo[, , tt]
        )
      }

      test_data <- subset(df_full, year == target_year)
      test_data <- test_data[
        complete.cases(test_data[, c(response, best_preds), drop = FALSE]),
      ]

      if (nrow(test_data) > 0) {
        preds_val <- predict(modelo_loyo, newdata = test_data)
        rmse <- sqrt(mean((test_data$f3 - preds_val)^2))
        denom_r2 <- sum((test_data$f3 - mean(test_data$f3))^2)
        r2 <- ifelse(
          denom_r2 > 0,
          1 - sum((test_data$f3 - preds_val)^2) / denom_r2,
          NA_real_
        )

        evaluation_log_loyo <- rbind(
          evaluation_log_loyo,
          data.frame(
            Month = mes,
            Model = paste("RF ventana - LOYO excluye", target_year),
            Formula = deparse(formula_modelo),
            RMSE = rmse,
            R2 = r2,
            stringsAsFactors = FALSE
          )
        )
      }
    }
    
    # ---------------------- Modelo común final ----------------------
    set.seed(SEED + mes)
    final_model_mes <- randomForest::randomForest(
      formula_modelo,
      data = df_full,
      ntree = 250,
      mtry = mtry_value,
      importance = TRUE
    )
    
    final_models[[mes]] <- final_model_mes
    final_predictors[[mes]] <- best_preds
    n_months_modelled <- n_months_modelled + 1L
    
    if (isTRUE(SAVE_MODELS)) {
      save(
        final_model_mes,
        best_preds,
        formula_modelo,
        file = file.path(output_dir_RData, paste0("RF_Model_", Modelo, "_", safe_biome_name, "_Mes", mes, ".RData"))
      )
    }
    
    # ---------------------- SHAP ----------------------
    if (length(best_preds) > 0 && !is.null(final_model_mes)) {
      
      set.seed(SEED + 100 + mes)
      X <- df_full[, best_preds, drop = FALSE]
      
      shap_values <- fastshap::explain(
        object = final_model_mes,
        X = X,
        pred_wrapper = predict,
        nsim = N_SHAP
      )
      
      shap_df <- as.data.frame(shap_values)
      
      shap_summary <- shap_df %>%
        dplyr::summarise(dplyr::across(dplyr::everything(), ~mean(abs(.), na.rm = TRUE))) %>%
        tidyr::pivot_longer(
          cols = dplyr::everything(),
          names_to = "Variable",
          values_to = "MeanAbsSHAP"
        ) %>%
        dplyr::mutate(Biome = safe_biome_name, Month = mes) %>%
        dplyr::arrange(dplyr::desc(MeanAbsSHAP)) %>%
        dplyr::mutate(Rank = dplyr::row_number())
      
      df_shap_global[[paste0(safe_biome_name, "_", mes)]] <- shap_summary
      
    } else {
      cat("[SHAP] No se puede calcular SHAP para mes", mes, "en", safe_biome_name, "\n")
      
      vars_na <- data.frame(
        Variable = best_preds,
        MeanAbsSHAP = NA_real_,
        Biome = safe_biome_name,
        Month = mes,
        Rank = NA_integer_
      )
      
      df_shap_global[[paste0(safe_biome_name, "_", mes)]] <- vars_na
    }
    
    # ---------------------- Predicción COMMON final ----------------------
    for (tt in mes_indices_central) {
      
      BA_FireHarmonized_common[, , tt] <- predict_layer_rf(
        modelo = final_model_mes,
        arr = arr_common,
        tt = tt,
        best_preds = best_preds,
        area_matrix_crop = area_matrix_crop,
        base_layer = BA_FireHarmonized_common[, , tt]
      )
    }
    
    preds_val <- predict(final_model_mes, newdata = df_full)
    rmse <- sqrt(mean((df_full$f3 - preds_val)^2))
    denom_r2 <- sum((df_full$f3 - mean(df_full$f3))^2)
    r2 <- ifelse(
      denom_r2 > 0,
      1 - sum((df_full$f3 - preds_val)^2) / denom_r2,
      NA_real_
    )
    
    evaluation_log <- rbind(
      evaluation_log,
      data.frame(
        Month = mes,
        Model = "RF ventana - todos los años COMMON",
        Formula = deparse(formula_modelo),
        RMSE = rmse,
        R2 = r2,
        stringsAsFactors = FALSE
      )
    )
    
    rm(df_full)
    gc()
  }
  
  # ---------------------------------------------------------------------------
  # 8.6 Guardar evaluación, predictores usados y SHAP
  # ---------------------------------------------------------------------------
  # Salidas LOYO comentadas por decisión metodológica actual.
  
  ev_loyo <- file.path(output_dir_csv, paste0("Evaluacion_RF_", safe_biome_name, "_LOYO.csv"))
  write.csv(evaluation_log_loyo, file = ev_loyo, row.names = FALSE)
  
  ev_common <- file.path(output_dir_csv, paste0("Evaluacion_RF_", safe_biome_name, "_COMMON.csv"))
  write.csv(evaluation_log, file = ev_common, row.names = FALSE)
  
  selected_used_csv <- file.path(output_dir_csv, paste0("SelectedPredictors_Used_", safe_biome_name, "_COMMON.csv"))
  write.csv(table_selected_used, file = selected_used_csv, row.names = FALSE)
  
  if (length(df_shap_global) > 0) {
    
    shap_all <- dplyr::bind_rows(df_shap_global, .id = "source")
    shap_csv <- file.path(output_dir_csv, paste0("SHAP_RF_", safe_biome_name, "_COMMON.csv"))
    write.csv(shap_all, file = shap_csv, row.names = FALSE)
    
    safe_biome_clean <- gsub("_", " ", gsub("__", " ", safe_biome_name))
    
    shap_grid <- tidyr::expand_grid(
      Variable = candidate_predictors,
      Month = 1:12
    ) %>%
      dplyr::left_join(
        shap_all %>% dplyr::select(Variable, Month, MeanAbsSHAP, Rank),
        by = c("Variable", "Month")
      ) %>%
      dplyr::mutate(
        Predictor = PREDICTOR_LABELS[Variable],
        MeanAbsSHAP_plot = MeanAbsSHAP,
        MeanAbsSHAP_plot = ifelse(MeanAbsSHAP_plot == 0, NA_real_, MeanAbsSHAP_plot)
      )
    
    pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(12)
    pal[1] <- "grey80"
    
    p_meanSHAP <- ggplot(
      shap_grid,
      aes(
        x = Month,
        y = factor(Predictor, levels = PREDICTOR_LABELS),
        fill = MeanAbsSHAP_plot
      )
    ) +
      geom_tile(color = "black", size = 0.2) +
      scale_x_continuous(
        breaks = 1:12,
        labels = month.abb,
        expand = c(0, 0)
      ) +
      scale_fill_stepsn(
        colours = pal,
        limits = c(0, 5),
        breaks = seq(0, 5, 0.5),
        oob = scales::squish,
        na.value = "grey80",
        name = "Mean\nAbs SHAP"
      ) +
      labs(x = "Month", y = "Predictor") +
      ggtitle(safe_biome_clean) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "white", fill = NA),
        legend.key.width = grid::unit(1.25, "cm"),
        legend.key.height = grid::unit(2, "cm"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 16)
      )
    
    jpeg(
      file = file.path(output_dir_plot_rle, paste0(safe_biome_name, "_meanSHAP.jpeg")),
      width = 2000,
      height = 1500,
      res = 300
    )
    print(p_meanSHAP)
    dev.off()
    
    max_rank_fixed <- length(candidate_predictors)
    rank_levels <- as.character(1:max_rank_fixed)
    pal_rank <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(max_rank_fixed)
    names(pal_rank) <- rank_levels
    
    rank_grid <- shap_grid %>%
      dplyr::mutate(
        Rank = factor(Rank, levels = rank_levels),
        Predictor = factor(Predictor, levels = PREDICTOR_LABELS)
      )
    
    dummy_levels <- tibble::tibble(
      Predictor = factor(PREDICTOR_LABELS[1], levels = PREDICTOR_LABELS),
      Month = 1,
      Rank = factor(rank_levels, levels = rank_levels)
    )
    
    rank_grid_full <- dplyr::bind_rows(rank_grid, dummy_levels)
    
    p_rank <- ggplot(
      rank_grid_full,
      aes(x = Month, y = Predictor, fill = Rank)
    ) +
      geom_tile(
        data = rank_grid,
        color = "black",
        size = 0.2
      ) +
      geom_tile(
        data = dummy_levels,
        alpha = 0,
        show.legend = TRUE
      ) +
      scale_x_continuous(
        breaks = 1:12,
        labels = month.abb,
        expand = c(0, 0)
      ) +
      scale_fill_manual(
        values = pal_rank,
        drop = FALSE,
        na.value = "grey80",
        name = "Rank"
      ) +
      labs(x = "Month", y = "Predictor") +
      ggtitle(safe_biome_clean) +
      theme_minimal() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "white", fill = NA),
        legend.key.width = grid::unit(1.0, "cm"),
        legend.key.height = grid::unit(1.0, "cm"),
        legend.title = element_text(size = 18, face = "bold"),
        legend.text = element_text(size = 14)
      )
    
    jpeg(
      file = file.path(output_dir_plot_rle, paste0(safe_biome_name, "_rankSHAP.jpeg")),
      width = 2000,
      height = 1500,
      res = 300
    )
    print(p_rank)
    dev.off()
  }
  
  # ---------------------------------------------------------------------------
  # 8.7 Armonización del periodo completo retrospectivo
  # ---------------------------------------------------------------------------
  
  for (mes in 1:12) {
    
    cat("\nArmonizando periodo completo para mes:", mes, "\n")
    
    if (is.null(final_models[[mes]]) || length(final_predictors[[mes]]) == 0) {
      cat("No hay modelo final para mes:", mes, ". Se conserva F51.\n")
      next
    }
    
    mes_indices_full <- which(lubridate::month(time_seq_full) == mes)
    
    for (tt in mes_indices_full) {
      
      BA_FireHarmonized_full[, , tt] <- predict_layer_rf(
        modelo = final_models[[mes]],
        arr = arr,
        tt = tt,
        best_preds = final_predictors[[mes]],
        area_matrix_crop = area_matrix_crop,
        base_layer = BA_FireHarmonized_full[, , tt]
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # 8.8 Series temporales, estadísticas y guardado por bioma
  # ---------------------------------------------------------------------------
  
  ba_total_Fire51 <- sapply(
    1:ntime_common,
    function(m) sum(obs_BA_Fire51_biome[, , m], na.rm = TRUE)
  )
  
  ba_total_FireS3 <- sapply(
    1:ntime_common,
    function(m) sum(obs_BA_FireS3_biome[, , m], na.rm = TRUE)
  )
  
  ba_Harmonized_loyo <- sapply(
    1:ntime_common,
    function(m) sum(BA_FireHarmonized_common_loyo[, , m], na.rm = TRUE)
  )
  
  ba_Harmonized_common <- sapply(
    1:ntime_common,
    function(m) sum(BA_FireHarmonized_common[, , m], na.rm = TRUE)
  )
  
  series_temporales_loyo <- data.frame(
    tiempo = 1:ntime_common,
    Date = time_seq_common,
    BA_FireCCI51 = ba_total_Fire51,
    BA_FireCCIS311 = ba_total_FireS3,
    BA_Harmonized = ba_Harmonized_loyo
  )
  
  archivo_series_loyo <- file.path(
    output_dir_csv,
    paste0("series_temporales_", Modelo, "_", safe_biome_name, "_LOYO.csv")
  )
  
#   write.csv(series_temporales_loyo, file = archivo_series_loyo, row.names = FALSE)
  
  series_temporales_common <- data.frame(
    tiempo = 1:ntime_common,
    Date = time_seq_common,
    BA_FireCCI51 = ba_total_Fire51,
    BA_FireCCIS311 = ba_total_FireS3,
    BA_Harmonized = ba_Harmonized_common
  )
  
  archivo_series_common <- file.path(
    output_dir_csv,
    paste0("series_temporales_", Modelo, "_", safe_biome_name, "_COMMON.csv")
  )
  
  write.csv(series_temporales_common, file = archivo_series_common, row.names = FALSE)
  
  stats_loyo <- rbind(
    cbind(Modelo = "F5_vs_S3", calc_stats(ba_total_FireS3, ba_total_Fire51)),
    cbind(Modelo = "Harmonized_LOYO_vs_S3", calc_stats(ba_total_FireS3, ba_Harmonized_loyo))
  )
  
  write.csv(
    stats_loyo,
    file = file.path(output_dir_csv, paste0("Statistics_temporal_series_", Modelo, "_", safe_biome_name, "_LOYO.csv")),
    row.names = FALSE
  )
  
  stats_common <- rbind(
    cbind(Modelo = "F5_vs_S3", calc_stats(ba_total_FireS3, ba_total_Fire51)),
    cbind(Modelo = "Harmonized_COMMON_vs_S3", calc_stats(ba_total_FireS3, ba_Harmonized_common))
  )
  
  write.csv(
    stats_common,
    file = file.path(output_dir_csv, paste0("Statistics_temporal_series_", Modelo, "_", safe_biome_name, "_COMMON.csv")),
    row.names = FALSE
  )
  
  plot_common_timeseries(
    bioma = bioma,
    safe_biome_name = safe_biome_name,
    dates_plot = time_seq_common,
    ba_ref = ba_total_FireS3,
    ba_f51 = ba_total_Fire51,
    ba_harm = ba_Harmonized_loyo,
    suffix = "LOYO",
    title_suffix = "(LOYO)"
  )
  
  plot_common_timeseries(
    bioma = bioma,
    safe_biome_name = safe_biome_name,
    dates_plot = time_seq_common,
    ba_ref = ba_total_FireS3,
    ba_f51 = ba_total_Fire51,
    ba_harm = ba_Harmonized_common,
    suffix = "COMMON",
    title_suffix = "(COMMON)",
    stats_text = format_stats_text(stats_common)
  )
  
  ba_Fire51_full <- sapply(
    seq_along(time_seq_full),
    function(m) sum(obs_BA_Fire51_biome_full[, , m], na.rm = TRUE)
  )
  
  ba_FireS3_full <- sapply(
    seq_along(time_seq_full),
    function(m) sum(obs_BA_FireS3_biome_full[, , m], na.rm = TRUE)
  )
  ba_FireS3_full[!(time_seq_full %in% time_seq_common)] <- NA_real_
  
  ba_Harmonized_full <- sapply(
    seq_along(time_seq_full),
    function(m) sum(BA_FireHarmonized_full[, , m], na.rm = TRUE)
  )
  
  series_temporales_full <- data.frame(
    tiempo = seq_along(time_seq_full),
    Date = time_seq_full,
    BA_FireCCI51 = ba_Fire51_full,
    BA_FireCCIS311 = ba_FireS3_full,
    BA_Harmonized_Full = ba_Harmonized_full
  )
  
  write.csv(
    series_temporales_full,
    file = file.path(output_dir_csv, paste0("series_temporales_", Modelo, "_", safe_biome_name, "_FULL.csv")),
    row.names = FALSE
  )
  
  if (isTRUE(SAVE_TIMESERIES_PLOTS)) {
    
    max_y_full <- max(ba_FireS3_full, ba_Fire51_full, ba_Harmonized_full, na.rm = TRUE) * 1.2
    if (!is.finite(max_y_full) || max_y_full <= 0) {
      max_y_full <- 1
    }
    
    jpeg(
      filename = file.path(output_dir_plot, paste0(safe_biome_name, "_Time_series_FULL.jpeg")),
      width = 2000,
      height = 1200,
      res = 300
    )
    
    par(mar = c(5, 5, 6, 2))
    plot(
      time_seq_full,
      ba_Fire51_full,
      type = "l",
      lwd = 2,
      col = "orange",
      xlab = "",
      ylab = expression("Burned Area (km"^2*")"),
      main = paste(bioma, "\nFireCCIS311, FireCCI51 and Harmonised FULL"),
      ylim = c(0, max_y_full),
      cex.lab = 1.2,
      cex.axis = 1.1,
      cex.main = 1.3,
      xaxt = "n"
    )
    
    lines(time_seq_full, ba_FireS3_full, col = "blue", lwd = 2)
    lines(time_seq_full, ba_Harmonized_full, col = "darkred", lwd = 2)
    axis.Date(
      1,
      at = seq(min(time_seq_full), max(time_seq_full), by = "12 months"),
      format = "%Y",
      las = 1,
      cex.axis = 0.8
    )
    grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
    legend(
      "topright",
      legend = c("FireCCIS311", "FireCCI51", "Harmonised"),
      col = c("blue", "orange", "darkred"),
      lwd = c(2, 2, 2),
      lty = c(1, 1, 1),
      ncol = 1,
      bty = "n",
      cex = 0.6,
      xpd = TRUE,
      seg.len = 1,
      inset = c(0, 0)
    )
    
    full_summary_text <- format_full_summary_text(
      dates_plot = time_seq_full,
      ba_f51 = ba_Fire51_full,
      ba_harm = ba_Harmonized_full,
      ba_s3 = ba_FireS3_full,
      stats_common = stats_common
    )
    
    if (!is.null(full_summary_text) && nzchar(full_summary_text)) {
      usr <- par("usr")
      text(
        x = usr[1],
        y = usr[4],
        labels = full_summary_text,
        adj = c(0, 1),
        cex = 0.42,
        col = "black",
        lineheight = 1.05
      )
    }
    
    dev.off()
  }
  
  obs_s3_rdata <- file.path(
    output_dir_RData,
    paste0("BA_", Modelo, "_", safe_biome_name, "_FireCCIS311.RData")
  )
  
  obs_f51_rdata <- file.path(
    output_dir_RData,
    paste0("BA_", Modelo, "_", safe_biome_name, "_FireCCI51.RData")
  )
  
  common_rdata <- file.path(
    output_dir_RData,
    paste0("BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Common.RData")
  )
  
  full_rdata <- file.path(
    output_dir_RData,
    paste0("BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Full.RData")
  )
  
  metadata_rdata <- file.path(
    output_dir_RData,
    paste0("BA_", Modelo, "_", safe_biome_name, "_Merge_Metadata.RData")
  )
  
  save(
    obs_BA_FireS3_biome,
    lon_vec_crop,
    lat_vec_crop,
    file = obs_s3_rdata
  )
  
  save(
    obs_BA_Fire51_biome,
    lon_vec_crop,
    lat_vec_crop,
    file = obs_f51_rdata
  )
  
  save(
    BA_FireHarmonized_common_loyo,
    lon_vec_crop,
    lat_vec_crop,
    file = file.path(output_dir_RData, paste0("BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Common_Loyo.RData"))
  )
  
  save(
    BA_FireHarmonized_common,
    lon_vec_crop,
    lat_vec_crop,
    file = common_rdata
  )
  
  save(
    BA_FireHarmonized_full,
    lon_vec_crop,
    lat_vec_crop,
    file = full_rdata
  )
  
  save(
    safe_biome_name,
    lon_idx,
    lat_idx,
    mask_lonlat,
    lon_vec_crop,
    lat_vec_crop,
    file = metadata_rdata
  )
  
  # ---------------------------------------------------------------------------
  # 8.9 Unión global diferida
  # ---------------------------------------------------------------------------
  # En modo paralelo no se escriben arrays globales desde los procesos hijos.
  # Cada bioma guarda sus RData locales y el proceso padre reconstruye los
  # mosaicos globales después de terminar todos los biomas.
  
  t_end <- Sys.time()
  
  cat("============================================================\n")
  cat("Bioma terminado:", bioma, "\n")
  cat("Meses modelizados:", n_months_modelled, "/ 12\n")
  cat("Fin:", as.character(t_end), "\n")
  cat("Duración segundos:", as.numeric(difftime(t_end, t_start, units = "secs")), "\n")
  cat("============================================================\n")
  
  rm(arr, arr_common)
  gc()
  
  list(
    Biome = safe_biome_name,
    Status = "OK",
    n_months_modelled = n_months_modelled,
    seconds = as.numeric(difftime(t_end, t_start, units = "secs")),
    selected_csv = selected_info$file,
    evaluation_loyo_csv = ev_loyo,
    evaluation_common_csv = ev_common,
    common_rdata = common_rdata,
    full_rdata = full_rdata,
    merge_metadata_rdata = metadata_rdata
  )
}

# ==============================================================================
# 9. WRAPPER SEGURO POR BIOMA
# ==============================================================================

safe_process_biome_harmonisation <- function(bioma) {
  
  safe_biome_name <- safe_name(as.character(bioma))
  
  tryCatch(
    {
      process_biome_harmonisation(bioma)
    },
    error = function(e) {
      
      error_file <- file.path(
        output_dir_log,
        paste0("ERROR_Harmonisation_", safe_biome_name, ".txt")
      )
      
      msg <- paste0(
        "ERROR EN BIOMA: ", bioma, "\n",
        "Fecha: ", as.character(Sys.time()), "\n",
        "Mensaje: ", conditionMessage(e), "\n"
      )
      
      writeLines(msg, error_file)
      
      list(
        Biome = safe_biome_name,
        Status = "ERROR",
        n_months_modelled = NA_integer_,
        seconds = NA_real_,
        selected_csv = NA_character_,
        evaluation_loyo_csv = NA_character_,
        evaluation_common_csv = NA_character_,
        common_rdata = NA_character_,
        full_rdata = NA_character_,
        merge_metadata_rdata = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
}

# ==============================================================================
# 10. EJECUCIÓN PARALELA POR BIOMA
# ==============================================================================

cat("============================================================\n")
cat("Inicio ejecución segunda parte desde SelectedPredictors CSV\n")
cat("Sistema operativo:", .Platform$OS.type, "\n")
cat("N_CORES solicitado:", N_CORES, "\n")
cat("Número de biomas:", length(biomas_unique), "\n")
cat("Directorio SelectedPredictors:", selected_predictors_dir2, "\n")
cat("============================================================\n")

t_global_start <- Sys.time()

if (.Platform$OS.type == "unix" && N_CORES > 1) {
  
  results_list <- parallel::mclapply(
    biomas_unique,
    safe_process_biome_harmonisation,
    mc.cores = N_CORES,
    mc.preschedule = FALSE
  )
  
} else {
  
  warning(
    "No se usa mclapply porque el sistema no es unix o N_CORES <= 1. ",
    "Ejecución secuencial."
  )
  
  results_list <- lapply(
    biomas_unique,
    safe_process_biome_harmonisation
  )
}

t_global_end <- Sys.time()

# ==============================================================================
# 11. RESUMEN FINAL Y GUARDADO GLOBAL
# ==============================================================================

results_df <- dplyr::bind_rows(
  lapply(results_list, as.data.frame)
)

results_df$global_start <- as.character(t_global_start)
results_df$global_end <- as.character(t_global_end)

status_csv <- file.path(
  output_dir_csv,
  "Harmonisation_From_SelectedPredictors_Run_Status.csv"
)

write.csv(
  results_df,
  file = status_csv,
  row.names = FALSE
)

print(results_df)

# ------------------------------------------------------------------------------
# 11.1 Reconstruir mosaicos globales de forma secuencial
# ------------------------------------------------------------------------------

global_BA_FireHarmonized_common_loyo <- array(
  NA_real_,
  dim = c(nrows, ncols, length(ind_common))
)

global_BA_FireHarmonized_common <- array(
  NA_real_,
  dim = c(nrows, ncols, length(ind_common))
)

global_BA_FireHarmonized_full <- array(
  NA_real_,
  dim = dim(BA_Fire51_tot)
)

merge_biome_to_global <- function(row_i) {
  
  status_i <- as.character(row_i$Status)
  
  if (!identical(status_i, "OK")) {
    return(invisible(FALSE))
  }
  
  meta_file <- as.character(row_i$merge_metadata_rdata)
  common_file <- as.character(row_i$common_rdata)
  full_file <- as.character(row_i$full_rdata)
  
  is_missing_file <- function(x) {
    length(x) == 0 || is.na(x[1]) || !nzchar(x[1]) || !file.exists(x[1])
  }
  
  if (
    is_missing_file(meta_file) ||
      is_missing_file(common_file) ||
      is_missing_file(full_file)
  ) {
    warning("Faltan archivos RData para unir el bioma: ", as.character(row_i$Biome))
    return(invisible(FALSE))
  }
  
  meta_env <- new.env(parent = emptyenv())
  load(meta_file, envir = meta_env)
  
  common_env <- new.env(parent = emptyenv())
  load(common_file, envir = common_env)
  
  full_env <- new.env(parent = emptyenv())
  load(full_file, envir = full_env)
  
  lon_idx_i <- meta_env$lon_idx
  lat_idx_i <- meta_env$lat_idx
  mask_i <- meta_env$mask_lonlat
  
  BA_common_i <- common_env$BA_FireHarmonized_common
  BA_full_i <- full_env$BA_FireHarmonized_full
  
  for (tt in seq_len(dim(BA_common_i)[3])) {
    tmp_global <- global_BA_FireHarmonized_common[lon_idx_i, lat_idx_i, tt]
    tmp_global[mask_i] <- BA_common_i[, , tt][mask_i]
    global_BA_FireHarmonized_common[lon_idx_i, lat_idx_i, tt] <<- tmp_global
  }
  
  for (tt in seq_len(dim(BA_full_i)[3])) {
    tmp_global <- global_BA_FireHarmonized_full[lon_idx_i, lat_idx_i, tt]
    tmp_global[mask_i] <- BA_full_i[, , tt][mask_i]
    global_BA_FireHarmonized_full[lon_idx_i, lat_idx_i, tt] <<- tmp_global
  }
  
  rm(meta_env, common_env, full_env, BA_common_i, BA_full_i)
  gc()
  invisible(TRUE)
}

cat("Uniendo resultados locales por bioma en mosaicos globales...\n")

for (ii in seq_len(nrow(results_df))) {
  merge_biome_to_global(results_df[ii, , drop = FALSE])
}

save(
  global_BA_FireHarmonized_common_loyo,
  lon_range,
  lat_range,
  file = file.path(output_dir_RData, paste0("BA_", Modelo, "_global_BA_FireHarmonized_Common_Loyo.RData"))
)

save(
  global_BA_FireHarmonized_common,
  lon_range,
  lat_range,
  file = file.path(output_dir_RData, paste0("BA_", Modelo, "_global_BA_FireHarmonized_Common.RData"))
)

save(
  global_BA_FireHarmonized_full,
  lon_range,
  lat_range,
  file = file.path(output_dir_RData, paste0("BA_", Modelo, "_global_BA_FireHarmonized_Full.RData"))
)

# ==============================================================================
# 12. SERIES Y ESTADÍSTICAS GLOBALES
# ==============================================================================

BA_Fire51_common <- BA_Fire51_tot[, , ind_common, drop = FALSE]
BA_FireS3_common <- BA_FireS3_tot[, , ind_common, drop = FALSE]

ntime_common <- length(ind_common)

ba_total_Fire51 <- sapply(
  1:ntime_common,
  function(m) sum(BA_Fire51_common[, , m], na.rm = TRUE)
)

ba_total_FireS3 <- sapply(
  1:ntime_common,
  function(m) sum(BA_FireS3_common[, , m], na.rm = TRUE)
)

ba_total_harmonized_loyo <- sapply(
  1:ntime_common,
  function(m) sum(global_BA_FireHarmonized_common_loyo[, , m], na.rm = TRUE)
)

ba_total_harmonized_common <- sapply(
  1:ntime_common,
  function(m) sum(global_BA_FireHarmonized_common[, , m], na.rm = TRUE)
)

series_global_common <- data.frame(
  tiempo = 1:ntime_common,
  Date = dates_common,
  BA_FireCCI51 = ba_total_Fire51,
  BA_FireCCIS311 = ba_total_FireS3,
#   BA_Harmonized_LOYO = ba_total_harmonized_loyo,
  BA_Harmonized_COMMON = ba_total_harmonized_common
)

write.csv(
  series_global_common,
  file = file.path(output_dir_csv, paste0("series_temporales_", Modelo, "_XGlobal_COMMON.csv")),
  row.names = FALSE
)

stats_global_common <- rbind(
  cbind(Modelo = "F5_vs_S3", calc_stats(ba_total_FireS3, ba_total_Fire51)),
#   cbind(Modelo = "Harmonized_LOYO_vs_S3", calc_stats(ba_total_FireS3, ba_total_harmonized_loyo)),
  cbind(Modelo = "Harmonized_COMMON_vs_S3", calc_stats(ba_total_FireS3, ba_total_harmonized_common))
)

write.csv(
  stats_global_common,
  file = file.path(output_dir_csv, paste0("Statistics_temporal_series_", Modelo, "_XGlobal_COMMON.csv")),
  row.names = FALSE
)

ba_total_Fire51_full <- sapply(
  seq_along(dates_full),
  function(m) sum(BA_Fire51_tot[, , m], na.rm = TRUE)
)

ba_total_FireS3_full <- sapply(
  seq_along(dates_full),
  function(m) sum(BA_FireS3_tot[, , m], na.rm = TRUE)
)
ba_total_FireS3_full[!(dates_full %in% dates_common)] <- NA_real_

ba_total_harmonized_full <- sapply(
  seq_along(dates_full),
  function(m) sum(global_BA_FireHarmonized_full[, , m], na.rm = TRUE)
)

series_global_full <- data.frame(
  tiempo = seq_along(dates_full),
  Date = dates_full,
  BA_FireCCI51 = ba_total_Fire51_full,
  BA_FireCCIS311 = ba_total_FireS3_full,
  BA_Harmonized_FULL = ba_total_harmonized_full
)

write.csv(
  series_global_full,
  file = file.path(output_dir_csv, paste0("series_temporales_", Modelo, "_XGlobal_FULL.csv")),
  row.names = FALSE
)

if (isTRUE(SAVE_TIMESERIES_PLOTS)) {
  
  max_y <- max(
    ba_total_FireS3,
    ba_total_Fire51,
    ba_total_harmonized_loyo,
    ba_total_harmonized_common,
    na.rm = TRUE
  ) * 1.2
  
  if (!is.finite(max_y) || max_y <= 0) {
    max_y <- 1
  }
  
  jpeg(
    filename = file.path(output_dir_plot, "XGLOBAL_COMMON_time_series.jpeg"),
    width = 2000,
    height = 1200,
    res = 300
  )
  
  par(mar = c(5, 5, 6, 2))
  plot(
    dates_common,
    ba_total_FireS3,
    type = "l",
    lwd = 2,
    col = "blue",
    xlab = "",
    ylab = expression("Burned Area (km"^2*")"),
    main = "Global COMMON: FireCCIS311, FireCCI51 and Harmonised",
    ylim = c(0, max_y),
    xaxt = "n"
  )
  lines(dates_common, ba_total_Fire51, col = "orange", lwd = 2)
#   lines(dates_common, ba_total_harmonized_loyo, col = "darkred", lwd = 2)
  lines(dates_common, ba_total_harmonized_common, col = "darkred", lwd = 2)
  axis.Date(
    1,
    at = seq(min(dates_common), max(dates_common), by = "3 months"),
    format = "%b-%Y",
    las = 1,
    cex.axis = 0.5
  )
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  legend(
    "topright",
    legend = c("FireCCIS311", "FireCCI51", "Harmonised COMMON"),
    col = c("blue", "orange", "darkred"),
    lwd = c(2, 2, 2),
    lty = c(1, 1, 1),
    bty = "n",
    cex = 0.6
  )
  
  global_common_stats_text <- format_stats_text(stats_global_common)
  
  if (!is.null(global_common_stats_text) && nzchar(global_common_stats_text)) {
    usr <- par("usr")
    text(
      x = usr[1],
      y = usr[4],
      labels = global_common_stats_text,
      adj = c(0, 1),
      cex = 0.42,
      col = "black",
      lineheight = 1.05
    )
  }
  
  dev.off()
  
  max_y_full <- max(ba_total_FireS3_full, ba_total_Fire51_full, ba_total_harmonized_full, na.rm = TRUE) * 1.2
  if (!is.finite(max_y_full) || max_y_full <= 0) {
    max_y_full <- 1
  }
  
  jpeg(
    filename = file.path(output_dir_plot, "XGLOBAL_FULL_time_series.jpeg"),
    width = 2000,
    height = 1200,
    res = 300
  )
  
  par(mar = c(5, 5, 6, 2))
  plot(
    dates_full,
    ba_total_Fire51_full,
    type = "l",
    lwd = 2,
    col = "orange",
    xlab = "",
    ylab = expression("Burned Area (km"^2*")"),
    main = "Global FULL: FireCCIS311, FireCCI51 and Harmonised",
    ylim = c(0, max_y_full),
    xaxt = "n"
  )
  lines(dates_full, ba_total_FireS3_full, col = "blue", lwd = 2)
  lines(dates_full, ba_total_harmonized_full, col = "darkred", lwd = 2)
  axis.Date(
    1,
    at = seq(min(dates_full), max(dates_full), by = "12 months"),
    format = "%Y",
    las = 1,
    cex.axis = 0.8
  )
  grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")
  legend(
    "topright",
    legend = c("FireCCIS311", "FireCCI51", "Harmonised FULL"),
    col = c("blue", "orange", "darkred"),
    lwd = c(2, 2, 2),
    lty = c(1, 1, 1),
    bty = "n",
    cex = 0.6
  )
  
  global_full_summary_text <- format_full_summary_text(
    dates_plot = dates_full,
    ba_f51 = ba_total_Fire51_full,
    ba_harm = ba_total_harmonized_full,
    ba_s3 = ba_total_FireS3_full,
    stats_common = stats_global_common
  )
  
  if (!is.null(global_full_summary_text) && nzchar(global_full_summary_text)) {
    usr <- par("usr")
    text(
      x = usr[1],
      y = usr[4],
      labels = global_full_summary_text,
      adj = c(0, 1),
      cex = 0.42,
      col = "black",
      lineheight = 1.05
    )
  }
  
  dev.off()
}

cat("============================================================\n")
cat("Ejecución terminada.\n")
cat("Inicio:", as.character(t_global_start), "\n")
cat("Fin:", as.character(t_global_end), "\n")
cat("Duración total:", as.numeric(difftime(t_global_end, t_global_start, units = "hours")), "horas\n")
cat("Resumen guardado en:\n")
cat(status_csv, "\n")
cat("Logs por bioma en:\n")
cat(output_dir_log, "\n")
cat("============================================================\n")

gc()
