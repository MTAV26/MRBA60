# ==============================================================================
# MRBA60 / FireCCI60
# SOLO selección de predictores por bioma y mes:
#   1) Filtro de autocorrelación
#   2) Recursive Feature Elimination
#
# Guarda:
#   SelectedPredictors_<bioma>_COMMON.csv
#
# NO ejecuta la armonización posterior.
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
  "caret",
  "parallel"
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
  library(caret)
  library(parallel)
})

sf::sf_use_s2(FALSE)

# ==============================================================================
# 1. CONFIGURACIÓN GENERAL
# ==============================================================================

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

output_dir     <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv <- file.path(output_dir, "csv")
output_dir_log <- file.path(output_dir, "logs_selected_predictors")

dirs <- c(output_dir, output_dir_csv, output_dir_log)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
}

# Número de núcleos: un bioma por núcleo
N_CORES <- min(3, max(1, parallel::detectCores() - 1))

SEED <- 123

cat("N_CORES =", N_CORES, "\n")

# Biomas a procesar.
# Usa NULL para todos, o vector numérico para solo algunos.
idx_biomas_pendientes <- unique(c(7, 8)) #NULL
# idx_biomas_pendientes <- unique(c(17:50, 9, 7, 8)) #NULL
# idx_biomas_pendientes <- NULL

# ==============================================================================
# 2. FUNCIONES AUXILIARES
# ==============================================================================
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_nc_var.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/safe_name.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/apply_mask_3d_2.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/build_prediction_df.R")



# ==============================================================================
# 3. CARGA DE DATOS GLOBALES
# ==============================================================================
load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))

dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

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
# 5. COORDENADAS
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


# 7. FUNCIÓN PRINCIPAL: SOLO SELECTED PREDICTORS
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/process_biome_selected_predictors.R")
# 8. WRAPPER SEGURO POR BIOMA
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/safe_process_biome_selected_predictors.R")

# ==============================================================================
# 9. EJECUCIÓN PARALELA POR BIOMA
# ==============================================================================

cat("============================================================\n")
cat("Inicio ejecución SOLO SelectedPredictors\n")
cat("Sistema operativo:", .Platform$OS.type, "\n")
cat("N_CORES solicitado:", N_CORES, "\n")
cat("Número de biomas:", length(biomas_unique), "\n")
cat("============================================================\n")

t_global_start <- Sys.time()

if (.Platform$OS.type == "unix" && N_CORES > 1) {
  
  results_list <- parallel::mclapply(
    biomas_unique,
    safe_process_biome_selected_predictors,
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
    safe_process_biome_selected_predictors
  )
}

t_global_end <- Sys.time()

# ==============================================================================
# 10. RESUMEN FINAL
# ==============================================================================

results_df <- dplyr::bind_rows(
  lapply(results_list, as.data.frame)
)

results_df$global_start <- as.character(t_global_start)
results_df$global_end <- as.character(t_global_end)

status_csv <- file.path(
  output_dir_csv,
  "SelectedPredictors_Run_Status.csv"
)

write.csv(
  results_df,
  file = status_csv,
  row.names = FALSE
)

print(results_df)

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
