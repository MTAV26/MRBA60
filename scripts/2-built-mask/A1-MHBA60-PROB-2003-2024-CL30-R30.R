# ======================================================================
# 0. LIMPIEZA Y PAQUETES
# ======================================================================
rm(list = ls())
graphics.off()
invisible(gc())

pkgs <- unique(c(
  "ggplot2", "ggtext", "ncdf4", "sp", "fields", "maps", "RColorBrewer",
  "dplyr", "lubridate", "MASS", "sf", "terra", "raster", "rworldmap",
  "graticule", "rnaturalearth", "rnaturalearthdata", "viridis", "caret",
  "gplots", "dendextend", "corrplot", "randomForest", "tidyr", "tibble",
  "cowplot", "fastshap", "pheatmap", "scales", "grid", "data.table", "stringr"
))
invisible(lapply(pkgs, function(x) library(x, character.only = TRUE)))

# ======================================================================
# 1. CONFIGURACIÓN GENERAL
# ======================================================================
Modelo <- "MRBA60-2003-2024-V1"

# Ruta base de datos de entrada
data_dir <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

# Producto de Active Fire a usar:
#   "FRPsum"  -> MODIS-FRPsum_conf30_angle30-200301-202412-025.RData
#   "FRPmean" -> MODIS-FRPmean_conf30_angle30-200301-202412-025.RData
active_fire_metric <- "AFcount"

# Directorios de salida
results_root <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask/"
output_dir <- file.path(results_root, Modelo)
# output_dir_csv <- file.path(output_dir, "csv/")
# output_dir_plot <- file.path(output_dir, "plot/")
# output_dir_plot_rle <- file.path(output_dir, "plot_rle/")
output_dir_RData <- file.path(output_dir, "RData/")

dirs <- c(output_dir,
          # output_dir_csv, output_dir_plot, output_dir_plot_rle,
          output_dir_RData)
for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ======================================================================
# 2. RUTAS DE ENTRADA
# ======================================================================
files_in <- list(
  longitude   = file.path(data_dir, "longitude.RData"),
  latitude    = file.path(data_dir, "latitude.RData"),
  firecci51   = file.path(data_dir, "FireCCI51_2003_2024_0.25degree.RData"),
  fireccis311 = file.path(data_dir, "FireCCIS311_2019_2024_0.25degree.RData"),
  active_fire = file.path(
    data_dir,
    paste0("MODIS-", active_fire_metric, "_conf30_angle30-200301-202412-025.RData")
  )
)

# Shapefile de biomas
file_biomes <- "/mnt/disco6tb/MRBA60/data/A1_RAW/MBC/continental-biomes_dinerstein_V10.shp"

# ======================================================================
# 3. COMPROBACIÓN DE ARCHIVOS
# ======================================================================
all_input_files <- c(unlist(files_in), file_biomes)
missing_files <- all_input_files[!file.exists(all_input_files)]

if (length(missing_files) > 0) {
  stop(
    "Faltan los siguientes archivos de entrada:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ======================================================================
# 4. HELPER PARA CARGAR .RData SIN ENSUCIAR EL ENTORNO
# ======================================================================
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_rdata.R")

# ======================================================================
# 5. CARGAR LONGITUD Y LATITUD
# ======================================================================
lon <- read_rdata(files_in$longitude)
lat <- read_rdata(files_in$latitude)

# ======================================================================
# 6. FECHAS: PERÍODO COMPLETO Y PERÍODO COMÚN
# ======================================================================
dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

anni <- 2003:2024
mesi <- rep(1:12, length(anni))
fechas <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
inicio <- which(fechas == as.Date("2003-01-01"))
fin    <- which(fechas == as.Date("2024-12-01"))

# ======================================================================
# 7. LEER BIOMAS
# ======================================================================
biomas_shp <- st_read(file_biomes, quiet = TRUE)
biomas_shp <- st_transform(biomas_shp, crs = 4326)

# ======================================================================
# 8. LEER DATOS DE ACTIVE FIRE, FIRECCIS311 Y FIRECCI51
# ======================================================================

# ----------------------------------------------------------------------
# 8.1 Active Fire
# NOTA:
# Se mantiene el nombre count_ActiveFire_tot por compatibilidad con el resto
# del script, aunque el archivo cargado pueda contener FRPsum o FRPmean.
# ----------------------------------------------------------------------
count_ActiveFire_tot <- read_rdata(files_in$active_fire)
count_ActiveFire_tot[count_ActiveFire_tot == 0] <- NA
invisible(gc())

# ----------------------------------------------------------------------
# 8.2 FireCCIS311 (2019-2024) -> respuesta f3
# ----------------------------------------------------------------------
BA_FireS3 <- read_rdata(files_in$fireccis311) / 1e6
BA_FireS3[BA_FireS3 == 0] <- NA
invisible(gc())

# ----------------------------------------------------------------------
# 8.3 FireCCI51 (2003-2024) -> predictor f5
# ----------------------------------------------------------------------
BA_Fire51_tot <- read_rdata(files_in$firecci51) / 1e6
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
invisible(gc())

nrows <- dim(BA_Fire51_tot)[1]
ncols <- dim(BA_Fire51_tot)[2]
time_common <- dim(BA_FireS3)[3]

# FireCCIS311 expandido a todo el periodo 2003-2024
BA_FireS3_tot <- array(NA_real_, dim = dim(BA_Fire51_tot))
BA_FireS3_tot[, , ind_common] <- BA_FireS3

# ======================================================================
# 10. MÁSCARA PARA EL ÁREA MÁXIMA DE CELDA (km²)
# ======================================================================
cell_area_constant <- (111.32 * 0.25)^2
area_by_row <- cell_area_constant * cos(lat * pi / 180)

nrow_grid <- length(lat)
ncol_grid <- length(lon)

area_matrix <- matrix(area_by_row, nrow = nrow_grid, ncol = ncol_grid, byrow = FALSE)
area_matrix <- t(area_matrix)

# ======================================================================
# 11. MÁSCARA PARA LIMITAR EL CÁLCULO A CELDAS CON INFORMACIÓN
# ======================================================================
# status_matrix2_tot = 1 si al menos una de las tres fuentes tiene valor > 0
# en esa celda/mes; en caso contrario, NA.

s_present  <- !is.na(count_ActiveFire_tot) & count_ActiveFire_tot > 0
f5_present <- !is.na(BA_Fire51_tot)        & BA_Fire51_tot > 0
f3_present <- !is.na(BA_FireS3_tot)        & BA_FireS3_tot > 0

status_matrix2_tot <- array(NA_real_, dim = dim(BA_Fire51_tot))
status_matrix2_tot[s_present | f5_present | f3_present] <- 1

rm(s_present, f5_present, f3_present)
invisible(gc())

# Si en las celdas con status = 1 hay NA, asignar 0
BA_Fire51_tot[status_matrix2_tot == 1 & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[status_matrix2_tot == 1 & is.na(BA_FireS3_tot)] <- 0
count_ActiveFire_tot[status_matrix2_tot == 1 & is.na(count_ActiveFire_tot)] <- 0

# ======================================================================
# 12. RECONSTRUIR LA MALLA DE COORDENADAS
# ======================================================================
ncol_raster <- length(lon)
nrow_raster <- length(lat)

lon_mat <- matrix(
  rep(lon, each = nrow_raster),
  nrow = nrow_raster,
  ncol = ncol_raster,
  byrow = FALSE
)

lat_mat <- matrix(
  rep(lat, times = ncol_raster),
  nrow = nrow_raster,
  ncol = ncol_raster,
  byrow = FALSE
)

lon_range <- lon_mat[1, ]
lat_range <- lat_mat[, 1]

# ======================================================================
# 13. CONFIGURACIÓN DE LOCALE PARA FECHAS
# ======================================================================
old_locale <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")

# ======================================================================
# 14. ARRAYS GLOBALES DE SALIDA
# ======================================================================
global_BA_FireHarmonized_common <- array(NA_real_, dim = dim(BA_FireS3))
global_BA_FireHarmonized_full   <- array(NA_real_, dim = dim(BA_Fire51_tot))

global_BA_FireHarmonized_prob_common <- array(NA_real_, dim = dim(BA_FireS3))
global_BA_FireHarmonized_prob_full   <- array(NA_real_, dim = dim(BA_Fire51_tot))

# ======================================================================
# 15. BIOMAS A EJECUTAR
# ======================================================================
biomas_unique <- unique(biomas_shp$cont_bm)
biomas_to_run <- biomas_unique

# ======================================================================
# 16. HELPERS
# ======================================================================
make_safe <- function(x) {
  x <- gsub(" ", "_", x)
  gsub("[^[:alnum:]_]", "", x)
}

month_window <- function(m) {
  c(((m + 10) %% 12) + 1, m, (m %% 12) + 1)
}

# ======================================================================
# 17. VARIABLES DEL MODELO
# ======================================================================
# Ya NO se leen predictores desde CSV. Se usan únicamente los datos cargados.
# Predictores permitidos: f5 (FireCCI51) y count_ActiveFire (Active Fire).
# Respuesta: f3 (FireS3).
valid_predictors <- c("f5", "count_ActiveFire")
response <- "f3"
dim(status_matrix2_tot)
image.plot(lon, lat, status_matrix2_tot[,,1])


for (bioma in biomas_to_run) {
  #===========================
  # Preparación por bioma
  #===========================
  df_shap_global <- list()
  cat("Procesando bioma:", bioma, "\n")
  safe_biome_name <- gsub(" ", "_", bioma)
  safe_biome_name <- gsub("[^[:alnum:]_]", "", safe_biome_name)
  
  # Selección del bioma y bounding box
  bioma_sel <- biomas_shp %>% filter(cont_bm == bioma)
  bbox <- st_bbox(bioma_sel)
  cat("Bounding box del bioma:", bbox, "\n")
  
  # Índices dentro del bbox
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    cat("No se encontraron celdas en la grilla para el bioma:", bioma, "\n")
    next
  }
  
  # ---- Recorte de datos ----
  BA_FireS3_crop_tot  <- BA_FireS3_tot[lon_idx,lat_idx,]
  BA_Fire51_crop_tot  <- BA_Fire51_tot[lon_idx,lat_idx,]
  count_ActiveFire_crop_tot <- count_ActiveFire_tot[lon_idx, lat_idx,]
  
  # Recortar matrices del período completo (totales)
  BA_Fire51_crop_full   <- BA_Fire51_tot[lon_idx, lat_idx, ]
  count_ActiveFire_crop_full <- count_ActiveFire_tot[lon_idx, lat_idx, ]
  
  # Extraer las coordenadas recortadas
  area_matrix_crop <- area_matrix[lon_idx, lat_idx]
  status_crop <- status_matrix2_tot[lon_idx, lat_idx, ]
  
  # Coordenadas y matrices auxiliares
  lon_mat_crop <- lon_mat[lat_idx, lon_idx]
  lat_mat_crop <- lat_mat[lat_idx, lon_idx]
  lon_vec_crop <- as.vector(lon_mat_crop)
  lat_vec_crop <- as.vector(lat_mat_crop)
  lon_vec <- sort(unique(as.vector(lon_mat_crop)))
  lat_vec <- sort(unique(as.vector(lat_mat_crop)))
  
  # ---- Crear objeto espacial ----
  grid_points_crop <- st_as_sf(data.frame(
    lon = as.vector(lon_mat_crop),
    lat = as.vector(lat_mat_crop)),
    coords = c("lon", "lat"), crs = 4326
  )
  
  # ---- Asignación de biomas ----
  inter <- st_intersects(grid_points_crop, biomas_shp)
  bioma_asignado <- sapply(inter, function(i) {
    if (length(i) == 0) return("Ninguno")
    return(biomas_shp$cont_bm[i[1]])
  })
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  # ---- Reasignar puntos sin bioma pero con al menos una presencia en el tiempo ----
  presencia_en_serie <- apply(status_crop, c(1,2), function(x) any(x == 1, na.rm = TRUE))
  grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie))
  
  puntos_sin_bioma_con_presencia <- grid_points_biomas %>%
    filter(bioma_final == "Ninguno" & presencia_serie)
  
  puntos_con_bioma <- grid_points_biomas %>%
    filter(bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    nearest_idx <- st_nearest_feature(puntos_sin_bioma_con_presencia, puntos_con_bioma)
    grid_points_biomas$bioma_final[
      which(grid_points_biomas$bioma_final == "Ninguno" & grid_points_biomas$presencia_serie)
    ] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  # ---- Crear máscara lógica para el bioma seleccionado ----
  mask <- grid_points_biomas$bioma_final == bioma
  mask_matrix <- matrix(mask, nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE)
  mask_matrix_sorted <- mask_matrix[order(lat_vec), order(lon_vec)]
  
  # ---- Aplicar máscara a toda la serie ----
  n_time <- dim(status_crop)[3]
  status_masked <- array(NA, dim = c(length(lon_vec), length(lat_vec), n_time))
  
  for (t in 1:n_time) {
    slice <- status_crop[,,t]
    slice_masked <- ifelse(t(mask_matrix_sorted), slice, NA)
    status_masked[,,t] <- slice_masked
  }
  
  mask_final <- apply(status_masked, c(1,2), mean, na.rm = TRUE)
  mask_final_clean <- mask_final==1
  mask_final_clean[is.na(mask_final_clean)]=FALSE
  
  dim_crop <- dim(BA_FireS3_crop_tot)
  
  # Construir arrays para el período completo (usar recortes de totales)
  ntime_full <- dim(BA_Fire51_crop_full)[3]
  obs_BA_Fire51_biome_full <- array(NA, dim = dim(BA_Fire51_crop_full))
  obs_count_ActiveFire_biome_full <- array(NA, dim = dim(count_ActiveFire_crop_full))
  
  dim_crop_tot <- dim(BA_FireS3_crop_tot)
  
  obs_BA_FireS3_biome_full <- array(NA, dim = dim_crop_tot)
  obs_BA_Fire51_biome_full <- array(NA, dim = dim_crop_tot)
  obs_count_ActiveFire_biome_full <- array(NA, dim = dim_crop_tot)
  
  for (t in 1:dim_crop_tot[3]) {
    slice_S3 <- BA_FireS3_crop_tot[,,t]
    slice_51 <- BA_Fire51_crop_tot[,,t]
    slice_af <- count_ActiveFire_crop_tot[,,t]
    
    slice_S3[!t(mask_matrix)] <- NA
    slice_51[!t(mask_matrix)] <- NA
    slice_af[!t(mask_matrix)] <- NA
    
    obs_BA_FireS3_biome_full[,,t] <- slice_S3
    obs_BA_Fire51_biome_full[,,t] <- slice_51
    obs_count_ActiveFire_biome_full[,,t] <- slice_af
  }
  
  rm(slice_S3, slice_51, slice_af)
  gc()
  
  obs_BA_FireS3_biome <- obs_BA_FireS3_biome_full[ , , ind_common]
  obs_BA_Fire51_biome <- obs_BA_Fire51_biome_full[ , , ind_common]
  obs_count_ActiveFire_biome <- obs_count_ActiveFire_biome_full[ , , ind_common]
  
  # Verificar que el bioma registra incendios en el período común (usamos S3)
  total_fire <- sum(obs_BA_FireS3_biome, na.rm = TRUE)
  if (total_fire == 0) {
    cat("El bioma", bioma, "no registra incendios en el periodo de escalado. Se guarda imagen con aviso.\n")
    jpeg(filename = file.path(output_dir_plot, paste0("time_series_", safe_biome_name, "_No_data.jpg")),
         width = 1600, height = 800, res = 150)
    plot(1, type="n", axes = FALSE, xlab = "", ylab = "", main = bioma)
    text(1, 1, "No se registran incendios para este bioma en el periodo de escalado.", cex = 1.5)
    dev.off()
    next
  }
  
  bin_FireS3_biome  <- obs_BA_FireS3_biome
  bin_FireS3_biome[bin_FireS3_biome>0]=1
  bin_Fire51_biome  <- obs_BA_Fire51_biome
  bin_Fire51_biome[bin_Fire51_biome>0]=1
  
  # ==============================================================================
  # 5. MODELOS POR MES (PERÍODO COMÚN) usando SOLO f5 y count_ActiveFire
  # ==============================================================================
  time_seq_common <- dates_common
  ntime_common    <- length(time_seq_common)
  
  # Arrays de salida (dominio recortado)
  dims <- dim(obs_BA_Fire51_biome)
  BA_FireHarmonized_common      <- array(NA, dim = dims)
  BA_FireHarmonized_prob_common <- array(NA, dim = dims)
  
  # === Predictores fijos por mes (sin CSV) ===
  fixed_preds_by_month <- setNames(replicate(12, c("f5","count_ActiveFire"), simplify = FALSE),
                                   as.character(1:12))
  
  # Listas para guardar modelos finales por mes
  final_models_ba   <- vector("list", 12)  # BA (regresión)
  final_models_prob <- vector("list", 12)  # Prob (clasificación)
  
  # Seguridad
  valid_predictors <- c("f5","count_ActiveFire")
  response <- "f3"
  
  for (mes in 1:12) {
    cat("\nProcesando mes (período común, ventana móvil):", mes, "\n")
    
    # Índices del mes central (para predecir/validar)
    mes_indices_central <- which(lubridate::month(time_seq_common) == mes)
    
    # Ventana móvil por índice (evita arrastrar dic-2022 cuando m==1)
    ntime_common <- length(time_seq_common)
    mes_indices_window <- sort(unique(
      c(mes_indices_central,
        mes_indices_central - 1,
        mes_indices_central + 1)
    ))
    mes_indices_window <- mes_indices_window[mes_indices_window >= 1 & mes_indices_window <= ntime_common]
    
    monthly_data <- list()
    
    # Paso 1: Recopilar datos históricos de TODA la ventana trimestral
    for (t in mes_indices_window) {
      f3_layer    <- obs_BA_FireS3_biome[,,t]
      f5_layer    <- obs_BA_Fire51_biome[,,t]
      count_layer <- obs_count_ActiveFire_biome[,,t]
      
      idx <- which(!is.na(f3_layer) & !is.na(f5_layer) & !is.na(count_layer))
      if (length(idx) > 0) {
        current_year <- lubridate::year(time_seq_common[t])
        key <- as.character(current_year)
        df_new <- data.frame(
          f3               = f3_layer[idx],
          f5               = f5_layer[idx],
          count_ActiveFire = count_layer[idx],
          year             = current_year
        )
        if (!is.null(monthly_data[[key]])) {
          monthly_data[[key]] <- rbind(monthly_data[[key]], df_new)
        } else {
          monthly_data[[key]] <- df_new
        }
      }
    }
    
    # Si no hay datos en la ventana, usar originales en el mes central
    if (length(monthly_data) == 0 || length(mes_indices_central) == 0) {
      cat("Mes", mes, ": Ventana sin datos. Usando valores originales para el mes central.\n")
      if (length(mes_indices_central) > 0) {
        BA_FireHarmonized_common[,,mes_indices_central]      <- obs_BA_Fire51_biome[,,mes_indices_central]
        BA_FireHarmonized_prob_common[,,mes_indices_central] <- bin_Fire51_biome[,,mes_indices_central]
      }
      next
    }
    
    df_full <- do.call(rbind, monthly_data)
    df_full <- df_full[complete.cases(df_full), ]
    
    # Verificar mínimo por año (>=30 celdas en ≥2 años)
    year_counts <- table(df_full$year)
    if (sum(year_counts >= 30) < 2) {
      cat("Mes", mes, ": Datos insuficientes en ventana. Usando valores originales para el mes central.\n")
      BA_FireHarmonized_common[,,mes_indices_central]      <- obs_BA_Fire51_biome[,,mes_indices_central]
      BA_FireHarmonized_prob_common[,,mes_indices_central] <- bin_Fire51_biome[,,mes_indices_central]
      next
    }
    
    # ======== PREDICTORES FIJOS (sin CSV) PARA ESTE MES ========
    rfe_preds_raw <- fixed_preds_by_month[[as.character(mes)]]
    rfe_preds <- intersect(valid_predictors, rfe_preds_raw)
    
    # Base mínima por robustez (forzamos 'f5')
    mandatory <- c("f5")
    preds_mes <- unique(c(mandatory, rfe_preds))
    preds_mes <- preds_mes[preds_mes %in% names(df_full)]
    
    if (length(preds_mes) == 0L) {
      warning("Mes ", mes, ": sin predictores válidos; intentando usar 'f5' si existe.")
      if ("f5" %in% names(df_full)) {
        preds_mes <- "f5"
      } else {
        preds_mes <- intersect(valid_predictors, names(df_full))
      }
    }
    
    # Subconjuntos de entrenamiento (BA y Prob usan los mismos predictores)
    cols_train <- unique(c(response, preds_mes, "year"))
    cols_train <- cols_train[cols_train %in% names(df_full)]
    df_full_sub <- df_full[, cols_train, drop = FALSE]
    
    # Comprobación de seguridad
    if (!(response %in% names(df_full_sub)) || length(setdiff(names(df_full_sub), c(response,"year"))) < 1) {
      cat("Mes", mes, ": no hay variables suficientes para entrenar. Uso originales (mes central).\n")
      BA_FireHarmonized_common[,,mes_indices_central]      <- obs_BA_Fire51_biome[,,mes_indices_central]
      BA_FireHarmonized_prob_common[,,mes_indices_central] <- bin_Fire51_biome[,,mes_indices_central]
      next
    }
    
    # Fórmulas y mtry
    rhs <- paste(setdiff(names(df_full_sub), c(response,"year")), collapse = " + ")
    formula_ba   <- as.formula(paste("f3 ~", rhs))
    # Probabilidad: binaria (1 si f3 > 0)
    df_full_sub$burn_bin <- factor(ifelse(df_full_sub$f3 > 0, 1, 0), levels = c(0,1))
    formula_prob <- as.formula(paste("burn_bin ~", rhs))
    
    # mtry_value <- max(2, floor(sqrt(length(setdiff(names(df_full_sub), c(response,"year"))))))
    mtry_value <- 1
    
    cat("Mes", mes, "- predictores usados:", paste(setdiff(names(df_full_sub), c(response,"year","burn_bin")), collapse = ", "),
        "| mtry =", mtry_value, "\n")
    
    # ENTRENAMIENTO (BA)
    model_ba <- randomForest(
      formula_ba, data = df_full_sub,
      ntree = 250, mtry = mtry_value, importance = TRUE
    )
    final_models_ba[[mes]] <- model_ba
    
    # ENTRENAMIENTO (Prob)
    model_prob <- randomForest(
      formula_prob, data = df_full_sub,
      ntree = 250, mtry = mtry_value, importance = TRUE
    )
    final_models_prob[[mes]] <- model_prob
    
    # PREDICCIÓN AL PERÍODO COMÚN (SOLO mes central) — BA + Prob
    preds_cols <- setdiff(names(df_full_sub), c(response, "year", "burn_bin"))
    for (t in mes_indices_central) {
      f5_layer    <- obs_BA_Fire51_biome[,,t]
      count_layer <- obs_count_ActiveFire_biome[,,t]
      idx <- which(!is.na(f5_layer) & !is.na(count_layer))
      if (length(idx) > 0) {
        pred_df <- data.frame(
          f5               = f5_layer[idx],
          count_ActiveFire = obs_count_ActiveFire_biome[,,t][idx]
        )
        pred_df <- pred_df[, preds_cols, drop = FALSE]
        
        # BA (regresión) — con límites físicos
        preds_ba <- predict(model_ba, newdata = pred_df)
        preds_ba <- pmax(preds_ba, 0)
        preds_ba <- pmin(preds_ba, area_matrix_crop[idx])
        
        layer_ba <- BA_FireHarmonized_common[,,t]
        if (all(is.na(layer_ba))) layer_ba <- array(NA, dim = dim(obs_BA_Fire51_biome[,,t]))
        layer_ba[idx] <- preds_ba
        BA_FireHarmonized_common[,,t] <- layer_ba
        
        # Prob (clasificación)
        preds_prob <- predict(model_prob, newdata = pred_df, type = "prob")[, "1"]
        layer_prob <- BA_FireHarmonized_prob_common[,,t]
        if (all(is.na(layer_prob))) layer_prob <- array(NA, dim = dim(obs_BA_Fire51_biome[,,t]))
        layer_prob[idx] <- preds_prob
        BA_FireHarmonized_prob_common[,,t] <- layer_prob
      } else {
        # Fallback directo si no hay pixeles válidos este mes
        BA_FireHarmonized_common[,,t]      <- obs_BA_Fire51_biome[,,t]
        BA_FireHarmonized_prob_common[,,t] <- bin_Fire51_biome[,,t]
      }
    }
    gc()
  }
  
  
  # ==============================================================================
  # 6. ARMONIZACIÓN DEL PERÍODO COMPLETO (FULL) CON LOS MODELOS FINALES
  # ==============================================================================
  time_seq_full <- dates_full
  ntime_full <- length(time_seq_full)
  BA_FireHarmonized_full       <- array(NA, dim = dim(obs_BA_Fire51_biome_full))
  BA_FireHarmonized_prob_full  <- array(NA, dim = dim(obs_BA_Fire51_biome_full))
  
  # Fallback binario para FULL (desde Fire51 full)
  bin_Fire51_biome_full <- obs_BA_Fire51_biome_full
  bin_Fire51_biome_full[bin_Fire51_biome_full > 0] <- 1
  bin_Fire51_biome_full[!is.na(bin_Fire51_biome_full) & bin_Fire51_biome_full <= 0] <- 0
  
  for (mes in 1:12) {
    cat("\nArmonizando período completo para mes:", mes, "\n")
    mes_indices_full <- which(lubridate::month(time_seq_full) == mes)
    
    # Si no hay modelos entrenados para este mes, fallback a Fire51 y su binario
    if (is.null(final_models_ba[[mes]]) || is.null(final_models_prob[[mes]])) {
      if (length(mes_indices_full) > 0) {
        BA_FireHarmonized_full[,,mes_indices_full]      <- obs_BA_Fire51_biome_full[,,mes_indices_full]
        BA_FireHarmonized_prob_full[,,mes_indices_full] <- bin_Fire51_biome_full[,,mes_indices_full]
      }
      next
    }
    
    # Columnas de predicción coherentes con las fórmulas entrenadas
    preds_cols_prob <- setdiff(all.vars(formula(final_models_prob[[mes]])), "burn_bin")
    preds_cols_ba   <- setdiff(all.vars(formula(final_models_ba[[mes]])),   "f3")
    
    for (t in mes_indices_full) {
      f5_layer    <- obs_BA_Fire51_biome_full[,,t]
      count_layer <- obs_count_ActiveFire_biome_full[,,t]
      idx <- which(!is.na(f5_layer) & !is.na(count_layer))
      if (length(idx) > 0) {
        pred_df <- data.frame(
          f5               = f5_layer[idx],
          count_ActiveFire = obs_count_ActiveFire_biome_full[,,t][idx]
        )
        
        # BA
        pred_df_ba <- pred_df[, preds_cols_ba, drop = FALSE]
        preds_ba <- predict(final_models_ba[[mes]], newdata = pred_df_ba)
        preds_ba <- pmax(preds_ba, 0)
        preds_ba <- pmin(preds_ba, area_matrix_crop[idx])
        layer_ba <- BA_FireHarmonized_full[,,t]
        if (all(is.na(layer_ba))) layer_ba <- array(NA, dim = dim(obs_BA_Fire51_biome_full[,,t]))
        layer_ba[idx] <- preds_ba
        BA_FireHarmonized_full[,,t] <- layer_ba
        
        # Prob
        pred_df_prob <- pred_df[, preds_cols_prob, drop = FALSE]
        preds_prob <- predict(final_models_prob[[mes]], newdata = pred_df_prob, type = "prob")[, "1"]
        layer_prob <- BA_FireHarmonized_prob_full[,,t]
        if (all(is.na(layer_prob))) layer_prob <- array(NA, dim = dim(obs_BA_Fire51_biome_full[,,t]))
        layer_prob[idx] <- preds_prob
        BA_FireHarmonized_prob_full[,,t] <- layer_prob
      } else {
        # Fallback si no hay pixeles válidos en este timestep
        BA_FireHarmonized_full[,,t]      <- f5_layer
        BA_FireHarmonized_prob_full[,,t] <- bin_Fire51_biome_full[,,t]
      }
    }
  }
  
  # ==============================================================================
  # 7. ENSAMBLAJE A MATRICES GLOBALES
  # ==============================================================================
  # Período común
  for (t in 1:dim(obs_BA_FireS3_biome)[3]) {
    biome_values <- BA_FireHarmonized_prob_common[,,t]    # Capa armonizada del bioma (común)
    tmp_global <- global_BA_FireHarmonized_prob_common[lon_idx, lat_idx, t]
    tmp_global[ mask_final_clean ] <- biome_values[ mask_final_clean ]
    global_BA_FireHarmonized_prob_common[lon_idx, lat_idx, t] <- tmp_global
  }
  # Período completo
  for (t in 1:dim(BA_Fire51_crop_full)[3]) {
    biome_values <- BA_FireHarmonized_prob_full[,,t]      # Capa armonizada del bioma (completo)
    tmp_global <- global_BA_FireHarmonized_prob_full[lon_idx, lat_idx, t]
    tmp_global[ mask_final_clean ] <- biome_values[ mask_final_clean ]
    global_BA_FireHarmonized_prob_full[lon_idx, lat_idx, t] <- tmp_global
  }
  
  # ==============================================================================
  # 8. GUARDADO POR BIOMA
  # ==============================================================================
  save(BA_FireHarmonized_prob_common, lon_vec_crop, lat_vec_crop, 
       file = paste0(output_dir_RData, "BA_PROB_", Modelo, "_", safe_biome_name, "_FireHarmonized_Common.RData"))
  save(BA_FireHarmonized_prob_full, lon_vec_crop, lat_vec_crop, 
       file = paste0(output_dir_RData, "BA_PROB_", Modelo, "_", safe_biome_name, "_FireHarmonized_Full.RData"))
}

# ==============================================================================
# 9. GUARDAR OBJETOS GLOBALES
# ==============================================================================
save(global_BA_FireHarmonized_prob_common, lon_range, lat_range,
     file = paste0(output_dir_RData, "BA_", Modelo, "global_BA_PROB_FireHarmonized_Common.RData"))
save(global_BA_FireHarmonized_prob_full, lon_range, lat_range,
     file = paste0(output_dir_RData, "BA_", Modelo, "global_BA_PROB_FireHarmonized_Full.RData"))


