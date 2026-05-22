
rm(list = ls())
graphics.off()
gc()

# ---------------------------
# Librerías
# ---------------------------
suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(lubridate)
  library(terra)
  library(raster)
  library(ggplot2)
  library(ggtext)
  library(ncdf4)
  library(sp)
  library(fields)
  library(maps)
  library(RColorBrewer)
  library(rworldmap)
  library(graticule)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(viridis)
  library(caret)
  library(gplots)
  library(dendextend)
  library(corrplot)
  library(randomForest)
  library(tidyr)
  library(tibble)
  library(cowplot)
  library(fastshap)
  library(qmap)
  library(lubridate)
})

# installed.packages("qmap")
# ---------------------------
# Configuración y rutas
# ---------------------------


Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
output_dir          <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)

output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
# output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_RData <- paste0(output_dir, "/RData/")
output_dir_plot_rle <- paste0(output_dir, "/plot_QM001M/")

# ---------------------------
# Cargar lon/lat y generar mallas
# ---------------------------
# load(file.path(dir_oss, "out-verifications-2019-2022_025/longitude.RData")) # objeto: longitude
# load(file.path(dir_oss, "out-verifications-2019-2022_025/latitude.RData"))  # objeto: latitude
load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))
lon_range <- as.vector(lon)
lat_range <- as.vector(lat)

# mallas: filas = lat, columnas = lon
lon_mat <- matrix(rep(lon_range, each = length(lat_range)), nrow = length(lat_range), ncol = length(lon_range), byrow = FALSE)
lat_mat <- matrix(rep(lat_range, times = length(lon_range)), nrow = length(lat_range), ncol = length(lon_range), byrow = FALSE)

# ---------------------------
# Fechas
# ---------------------------
dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)  # 2019–2022
months_vec   <- month(dates_full)

# # ---------------------------
# # Cargar FireCCIS311 (solo periodo común 2019-2022) y normalizar unidades
# # ---------------------------
# load(paste0(dir_oss, "out-verifications-2019-2022_025/FireCCIS311_2019_2022_0.25degree-download.RData"))
# BA_FireS3 <- FireS3 / 1e6
# BA_FireS3[BA_FireS3 == 0] <- NA
# rm(FireS3); gc()
# # image(BA_FireS3[,,1])
# # ---------------------------
# # Cargar FireCCI51 total (2001-2022), recortar a 2003-2022 y normalizar
# # ---------------------------
# load(paste0(dir_oss, "out-verifications-2019-2022_025/FireCCI51_2001_2022_0.25degree-download.RData"))
# BA_Fire51_tot <- Fire51 / 1e6
# BA_Fire51_tot <- BA_Fire51_tot[,,25:264]  # 2003-01 a 2022-12
# BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
# rm(Fire51); gc()
# image(BA_Fire51_tot[,,1])



load(file.path(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))
BA_FireS3 <- s3 / 1e6
BA_FireS3[BA_FireS3 == 0] <- NA
rm(s3)
gc()

load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))
BA_Fire51_tot <- f51 / 1e6
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
rm(f51)
gc()


# load("/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/M3_RF-PROB-2003-2022-CL60-R30/RData/BA_M3_RF-PROB-2003-2022-CL60-R30global_BA_FireHarmonized_Full_Filtered_Restored.RData")
load(file.path(output_dir_RData,"BA_MRBA60.RData"))
# # save(BA_qm_global, file = outfile)
# dim(BA_qm_global_co)
# dim(BA_FireS3)
# BA_FIRE60
BA_FIRE60[BA_FIRE60 == 0] <- NA
BA_harmonised=BA_FIRE60
rm(BA_FIRE60)
gc()
image(BA_harmonised[,,1])
# BA_harmonised[is.na(BA_harmonised)]=0
# BA_final[is.na(BA_final)]=0
# BA_Fire51_tot=BA_final

nrows <- dim(BA_harmonised)[1]
ncols <- dim(BA_harmonised)[2]
time_common <- dim(BA_FireS3)[3]  # número de meses para el período común

BA_FireS3_tot <- array(NA, dim = c(nrows, ncols, dim(BA_harmonised)[3]))
BA_FireS3_tot[,,ind_common] <- BA_FireS3
# BA_FireS3_tot[BA_FireS3_tot==0]=NA



# image.plot(lon, lat, BA_Fire51_tot[,,1])

# ============================================================================
# MÁSCARA PARA EL ÁREA MÁXIMA DE CELDA (km²)
# ============================================================================
# cell_area_constant <- (111.32 * 0.25)^2  # en el ecuador
# area_by_row <- cell_area_constant * cos(lat * pi/180)

cell_area_constant <- (110.57 * 0.25) * (111.32 * 0.25)  # ≈ 769.29 km² en el ecuador
area_by_row <- cell_area_constant * cos(lat * pi/180)

nrow_grid <- length(lat)
ncol_grid <- 1440
area_matrix <- matrix(area_by_row, ncol = ncol_grid, nrow = nrow_grid, byrow = FALSE)
area_matrix <- t(area_matrix)
image(lon, lat, area_matrix)


nlon <- length(lon)
nlat <- length(lat)
# Matriz 2D de longitudes (constante respecto a lat): [lon x lat]
lon_matrix <- matrix(lon, nrow = nlon, ncol = nlat, byrow = FALSE)
image(lon, lat, lon_matrix, xlab = "lon", ylab = "lat",
      main = "Longitude (°) — matrix [lon × lat]")


# load("/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30/MASK_FIRE_PREPROC/status_matrix2_tot_filter.RData")



# ============================================================================
# MÁSCARA PARA LIMITAR EL CÁLCULO A CELDAS CON INFORMACIÓN
# ============================================================================
status_matrix2_tot <- array(NA, dim = c(nrows, ncols, dim(BA_Fire51_tot)[3]))
for (t in 1:dim(BA_Fire51_tot)[3]) {
  for (i in 1:nrows) {
    for (j in 1:ncols) {
      s_punt <- BA_harmonised[i, j, t]
      f_val  <- BA_Fire51_tot[i, j, t]
      f_val3 <- BA_FireS3_tot[i, j, t]
      
      s_present_punt <- (!is.na(s_punt)) && (s_punt > 0)
      f_present      <- (!is.na(f_val)) && (f_val > 0)
      f_present3     <- (!is.na(f_val3)) && (f_val3 > 0)
      
      if (s_present_punt && f_present && f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (s_present_punt && f_present && !f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (!s_present_punt && f_present && f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (s_present_punt && !f_present && f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (!s_present_punt && f_present && !f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (!s_present_punt && !f_present && f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else if (s_present_punt && !f_present && !f_present3) {
        status_matrix2_tot[i, j, t] <- 1
      } else {
        status_matrix2_tot[i, j, t] <- NA
      }
    }
  }
}
rm(f_present3, f_present, s_present_punt, f_val3, f_val, s_punt)
gc()

BA_Fire51_tot[status_matrix2_tot == 1 & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[status_matrix2_tot == 1 & is.na(BA_FireS3_tot)] <- 0
BA_harmonised[status_matrix2_tot == 1 & is.na(BA_harmonised)] <- 0

rm(BA_Fire51_tot)
gc()

image(status_matrix2_tot[,,1])
# BA_final[is.na(BA_final)]=0
BA_Fire51_tot=BA_harmonised

global_rmse_abs_full <- array(NA, dim = dim(BA_Fire51_tot))
global_rmse_rel_full    <- array(NA, dim = dim(BA_Fire51_tot))

# biomas_shp <- st_read("/mnt/disco6tb/FireCCI60/data_025/continental-biomes_dinerstein_V9/continental-biomes_dinerstein_V10.shp")

biomas_shp <- st_read(
  file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"),
  quiet = TRUE
)

biomas_shp <- st_transform(biomas_shp, crs = 4326)
biomas_unique <- unique(biomas_shp$cont_bm)
# biomas_shp <- st_read(
#   file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"),
#   quiet = TRUE
# 
# )



nbins    <- 10
deg_poly <- 1


# for (bioma in biomas_unique[1]) {  # para probar uno
for (bioma in biomas_unique) {
  df_shap_global <- list()
  cat("Procesando bioma:", bioma, "\n")
  safe_biome_name <- gsub(" ", "_", bioma)
  safe_biome_name <- gsub("[^[:alnum:]_]", "", safe_biome_name)
  
  # ---- Selección del bioma y bounding box ----
  bioma_sel <- biomas_shp %>% dplyr::filter(cont_bm == bioma)
  bbox <- sf::st_bbox(bioma_sel)
  cat("Bounding box del bioma:", bbox, "\n")
  
  # ---- Índices dentro del bbox ----
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    cat("No se encontraron celdas en la grilla para el bioma:", bioma, "\n")
    next
  }
  
  # ---- Recorte de datos ----
  BA_FireS3_crop_tot <- BA_FireS3_tot[lon_idx, lat_idx, ]
  BA_Fire51_crop_tot <- BA_Fire51_tot[lon_idx, lat_idx, ]
  
  BA_Fire51_crop_full <- BA_Fire51_tot[lon_idx, lat_idx, ]
  area_matrix_crop    <- area_matrix[lon_idx, lat_idx]
  lon_matrix_crop     <- lon_matrix[lon_idx, lat_idx]
  status_crop         <- status_matrix2_tot[lon_idx, lat_idx, ]
  
  # ---- Coordenadas y matrices auxiliares ----
  lon_mat_crop <- lon_mat[lat_idx, lon_idx]  # filas = lat, cols = lon
  lat_mat_crop <- lat_mat[lat_idx, lon_idx]
  lon_vec <- sort(unique(as.vector(lon_mat_crop)))
  lat_vec <- sort(unique(as.vector(lat_mat_crop)))
  
  # ---- Crear objeto espacial ----
  grid_points_crop <- sf::st_as_sf(
    data.frame(lon = as.vector(lon_mat_crop), lat = as.vector(lat_mat_crop)),
    coords = c("lon", "lat"), crs = 4326
  )
  
  # ---- Asignación de biomas ----
  inter <- sf::st_intersects(grid_points_crop, biomas_shp)
  bioma_asignado <- sapply(inter, function(i) if (length(i) == 0) "Ninguno" else biomas_shp$cont_bm[i[1]])
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  # ---- Reasignar puntos sin bioma pero con presencia temporal ----
  presencia_en_serie <- apply(status_crop, c(1, 2), function(x) any(x == 1, na.rm = TRUE))
  grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie))
  
  puntos_sin_bioma_con_presencia <- dplyr::filter(grid_points_biomas, bioma_final == "Ninguno" & presencia_serie)
  puntos_con_bioma               <- dplyr::filter(grid_points_biomas, bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    nearest_idx <- sf::st_nearest_feature(puntos_sin_bioma_con_presencia, puntos_con_bioma)
    grid_points_biomas$bioma_final[
      which(grid_points_biomas$bioma_final == "Ninguno" & grid_points_biomas$presencia_serie)
    ] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  # ---- Máscara lógica para el bioma seleccionado ----
  mask <- grid_points_biomas$bioma_final == bioma
  # mask es vectorizado en orden (lat, lon); para aplicarlo sobre [lon,lat] usamos t(mask_matrix)
  mask_matrix <- matrix(mask, nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE)
  
  # ---- Aplicar máscara a toda la serie status (si lo necesitas para depuración) ----
  n_time <- dim(status_crop)[3]
  status_masked <- array(NA, dim = c(length(lon_vec), length(lat_vec), n_time))
  for (t in 1:n_time) {
    slice <- status_crop[, , t]
    slice_masked <- ifelse(t(mask_matrix), slice, NA)
    status_masked[, , t] <- slice_masked
  }
  
  # ---- Construir arrays para el período completo con máscara ----
  dim_crop_tot <- dim(BA_FireS3_crop_tot)
  obs_BA_FireS3_biome_full <- array(NA_real_, dim = dim_crop_tot)  # [lon, lat, time] (mismo que crop_tot)
  obs_BA_Fire51_biome_full <- array(NA_real_, dim = dim_crop_tot)
  
  for (t in 1:dim_crop_tot[3]) {
    slice_S3 <- BA_FireS3_crop_tot[, , t]
    slice_51 <- BA_Fire51_crop_tot[, , t]
    slice_S3[!t(mask_matrix)] <- NA
    slice_51[!t(mask_matrix)] <- NA
    obs_BA_FireS3_biome_full[, , t] <- slice_S3
    obs_BA_Fire51_biome_full[, , t] <- slice_51
  }
  rm(slice_S3, slice_51); gc()
  
  obs_BA_FireS3_biome <- obs_BA_FireS3_biome_full[, , ind_common]
  obs_BA_Fire51_biome <- obs_BA_Fire51_biome_full[, , ind_common]
  
  # Producto predictor para el ajuste (tu armonizado a partir de Fire51 en el bioma)
  pred_harm_biome_full <- obs_BA_Fire51_biome_full
  
  # ---- Chequeo de incendios en el periodo común ----
  total_fire <- sum(obs_BA_FireS3_biome, na.rm = TRUE)
  if (total_fire == 0) {
    cat("El bioma", bioma, "no registra incendios en el periodo de escalado. Se guarda imagen con aviso.\n")
    jpeg(filename = file.path(output_dir_plot, paste0("time_series_", safe_biome_name, "_No_data.jpg")),
         width = 1600, height = 800, res = 150)
    plot(1, type = "n", axes = FALSE, xlab = "", ylab = "", main = bioma)
    text(1, 1, "No se registran incendios para este bioma en el periodo de escalado.", cex = 1.5)
    dev.off()
    next
  }
  
  # ==== Dimensiones y preasignación (ANTES de usar rmse_*_full) ====
  nlon       <- dim(obs_BA_FireS3_biome_full)[1]
  nlat       <- dim(obs_BA_FireS3_biome_full)[2]
  ntime_full <- dim(obs_BA_FireS3_biome_full)[3]
  
  rmse_abs_full <- array(NA_real_, dim = c(nlon, nlat, ntime_full))
  rmse_rel_full <- array(NA_real_, dim = c(nlon, nlat, ntime_full))
  
  # ==== Índices de entrenamiento y de predicción ====
  train_idx <- ind_common
  full_idx  <- seq_len(ntime_full)
  pre_idx   <- setdiff(full_idx, train_idx)
  
  # ==== Vectorizar SOLO el periodo común (train) ====
  pred_common_vec <- unlist(lapply(train_idx, function(t) as.vector(pred_harm_biome_full[, , t])))
  obs_common_vec  <- unlist(lapply(train_idx, function(t) as.vector(obs_BA_FireS3_biome_full[, , t])))
  
  ok <- is.finite(pred_common_vec) & is.finite(obs_common_vec)
  pred_common_vec <- pred_common_vec[ok]
  obs_common_vec  <- obs_common_vec[ok]
  
  # ==== Bins por cuantiles con SOLO train ====
  probs <- seq(0, 1, length.out = nbins + 1)
  brks  <- unique(quantile(pred_common_vec, probs = probs, na.rm = TRUE, type = 7))
  if (length(brks) < 2) {
    r <- range(pred_common_vec, na.rm = TRUE)
    if (r[1] == r[2]) r[2] <- r[2] + .Machine$double.eps
    brks <- r
  }
  brks <- sort(brks)
  
  bin_id <- cut(pred_common_vec, breaks = brks, include.lowest = TRUE, labels = FALSE, right = TRUE)
  
  df_bin <- data.frame(pred = pred_common_vec, obs = obs_common_vec, bin = bin_id) |>
    dplyr::group_by(bin) |>
    dplyr::summarise(
      pred_mean = mean(pred, na.rm = TRUE),
      rmse_abs  = sqrt(mean((pred - obs)^2, na.rm = TRUE)),
      obs_mean  = mean(obs, na.rm = TRUE),
      n         = dplyr::n(),
      .groups   = "drop"
    ) |>
    dplyr::filter(!is.na(bin))
  
  df_bin$rmse_rel <- 100 * df_bin$rmse_abs / pmax(df_bin$obs_mean, .Machine$double.eps)
  
  print(df_bin)
  # ==== Ajuste SOLO con train ====
  fit_abs <- lm(rmse_abs ~ pred_mean, data = df_bin)   # Puedes cambiar por Gamma/log si lo prefieres
  fit_rel <- lm(rmse_rel ~ pred_mean, data = df_bin)
  
  # Rango de pred_mean visto en TRAIN (para clip al predecir)
  # rng <- range(df_bin$pred_mean, na.rm = TRUE)
  
  # ==== Predicción en el periodo común ====
  for (t_full in train_idx) {
    v_pred <- as.vector(pred_harm_biome_full[, , t_full])
    # v_pred <- pmin(pmax(v_pred, rng[1]), rng[2])   # evitar extrapolación rara
    
    y_abs <- pmax(as.numeric(predict(fit_abs, newdata = data.frame(pred_mean = v_pred))), 0)
    y_rel <- pmax(as.numeric(predict(fit_rel, newdata = data.frame(pred_mean = v_pred))), 0)
    
    
    rmse_abs_full[, , t_full] <- matrix(y_abs, nrow = nlon, ncol = nlat)
    rmse_rel_full[, , t_full] <- matrix(y_rel, nrow = nlon, ncol = nlat)
  }
  
  # ==== Propagación al resto (2003–2018 y años fuera de train) ====
  for (t in pre_idx) {
    v_pred <- as.vector(pred_harm_biome_full[, , t])
    # v_pred <- pmin(pmax(v_pred, rng[1]), rng[2])
    
    y_abs <- as.numeric(predict(fit_abs, newdata = data.frame(pred_mean = v_pred)))
    y_rel <- as.numeric(predict(fit_rel, newdata = data.frame(pred_mean = v_pred)))
    
    rmse_abs_full[, , t] <- matrix(y_abs, nrow = nlon, ncol = nlat)
    rmse_rel_full[, , t] <- matrix(y_rel, nrow = nlon, ncol = nlat)
  }
  
  # ==== Volcado a cubos globales usando la MISMA máscara ====
  for (t in 1:dim(BA_Fire51_crop_full)[3]) {
    tmp_global <- global_rmse_abs_full[lon_idx, lat_idx, t]
    tmp_global[t(mask_matrix)] <- rmse_abs_full[, , t][t(mask_matrix)]
    global_rmse_abs_full[lon_idx, lat_idx, t] <- tmp_global
  }
  for (t in 1:dim(BA_Fire51_crop_full)[3]) {
    tmp_global <- global_rmse_rel_full[lon_idx, lat_idx, t]
    tmp_global[t(mask_matrix)] <- rmse_rel_full[, , t][t(mask_matrix)]
    global_rmse_rel_full[lon_idx, lat_idx, t] <- tmp_global
  }
  
  # (Opcional) quicklook
  image.plot(lon, lat, global_rmse_abs_full[,,13])
  print(summary(as.vector(global_rmse_abs_full)))
}


global_rmse_rel_full[global_rmse_rel_full<=0]=0
image.plot(lon, lat, global_rmse_rel_full[,,8])
print(summary(as.vector(global_rmse_rel_full)))

# global_rmse_rel_full[global_rmse_abs_full<=0]=0
image.plot(lon, lat, global_rmse_abs_full[,,8])
print(summary(as.vector(global_rmse_rel_full)))
# load("/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-2/RData/")

outfile <- file.path(output_dir_RData, "BA_Incertidumbre_HARMONISED_abs.RData")
save(global_rmse_abs_full, file = outfile)

outfile <- file.path(output_dir_RData, "BA_Incertidumbre_HARMONISED_rel.RData")
save(global_rmse_rel_full, file = outfile)


