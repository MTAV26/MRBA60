# ===============================
# CONFIGURACIÓN PREVIA
# ===============================
rm(list = ls())
graphics.off()
gc()

# Cargar librerías
library(ggplot2)
library(ggtext)
library(ncdf4)
library(sp)
library(fields)
library(maps)
library(RColorBrewer)
library(dplyr)
library(lubridate)
library(MASS)
library(sf)
library(terra)
library(raster)
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
library(ggplot2)
library(cowplot)
library(fastshap)
Modelo <- "B1-MRBA60-2003-2024"

# Directorios de trabajo y salida
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
# /mnt/disco6tb/MRBA60-2/results/

output_dir <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_plot_sc <- paste0(output_dir, "/plot_scatter/")
output_dir_RData <- paste0(output_dir, "/RData/")

# Crear directorios si no existen
dirs <- c(output_dir, output_dir_csv, output_dir_plot, output_dir_plot_rle, output_dir_RData, output_dir_plot_sc)
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# ============================================================================
# CARGAR LONGITUD Y LATITUD
# ============================================================================
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")

# ============================================================================
# FECHAS: PERÍODO COMPLETO Y PERÍODO COMÚN
# ============================================================================
dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

anni <- 2000:2024
mesi <- rep(1:12, length(anni))
fechas <- seq(as.Date("2000-01-01"), as.Date("2024-12-01"), by = "month")
inicio <- which(fechas == as.Date("2003-01-01"))
fin <- which(fechas == as.Date("2024-12-01"))

# ============================================================================
# LEER BIOMAS
# ============================================================================

biomas_shp <- st_read(file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"))
biomas_shp <- st_transform(biomas_shp, crs = 4326)
# # ============================================================================
# # LEER DATOS DE ACTIVE FIRE, FIRECCIS311 Y FIRECCI51
# # ============================================================================
# # Datos de Active Fire (se usan para el período común y completo)
# load(paste0(dir_oss, "MODIS-count_AvtiveFire_conf70-2001-2024-025_conf70.RData"))
# count_ActiveFire_tot <- count_AvtiveFire_conf70[,,25:264]
# # count_ActiveFire <- count_ActiveFire_tot[,,ind_common]
# count_ActiveFire_tot[count_ActiveFire_tot == 0] <- NA
# # count_ActiveFire[count_ActiveFire == 0] <- NA
# 
# FireCCIS311 (disponible solo para el período común)
load(paste0(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))
BA_FireS3 <- s3 / 1e6
# summary(as.vector(BA_FireS3))
BA_FireS3[BA_FireS3 == 0] <- NA
rm(s3)
gc()

# FireCCI51: se carga el total y se extrae el período común
load(paste0(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))
BA_Fire51_tot <- f51 / 1e6
BA_Fire51_tot <- BA_Fire51_tot#[,,25:264]
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
rm(f51)
gc()
# 
# 
nrows <- length(lon)
ncols <-length(lat)
time_common <- 264  # número de meses para el período común
# 
# # Crear arreglo para BA_FireS3_tot basado en el período común (para mantener la misma estructura)
BA_FireS3_tot <- array(NA, dim = c(nrows, ncols, dim(BA_Fire51_tot)[3]))
BA_FireS3_tot[,,ind_common] <- BA_FireS3

# ============================================================================
# MÁSCARA PARA EL ÁREA MÁXIMA DE CELDA (km²)
# ============================================================================
cell_area_constant <- (111.32 * 0.25)^2  # en el ecuador
area_by_row <- cell_area_constant * cos(lat * pi/180)
nrow_grid <- length(lat)
ncol_grid <- 1440
area_matrix <- matrix(area_by_row, ncol = ncol_grid, nrow = nrow_grid, byrow = FALSE)
area_matrix <- t(area_matrix)

# image(area_matrix)
load(paste0(dir_oss,"FireMask_AF3030F.RData"))
#dim(FireMask_AF3030F)
#image.plot(lon, lat, FireMask_AF3030F[,,2])
status_matrix2_tot=FireMask_AF3030F
rm(FireMask_AF3030F)
gc()
# # Si en las celdas de BA_Fire51 hay NA y status_matrix2 es 1, asignar 0
BA_Fire51_tot[status_matrix2_tot == 1 & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[status_matrix2_tot == 1 & is.na(BA_FireS3_tot)] <- 0
# count_ActiveFire_tot[status_matrix2_tot == 1 & is.na(count_ActiveFire_tot)] <- 0
# ============================================================================
# RECONSTRUIR LA MALLA DE COORDENADAS
# ============================================================================
ncol_raster <- length(lon)
nrow_raster <- length(lat)
lon_mat <- matrix(rep(lon, each = nrow_raster), nrow = nrow_raster, ncol = ncol_raster, byrow = FALSE)
lat_mat <- matrix(rep(lat, times = ncol_raster), nrow = nrow_raster, ncol = ncol_raster, byrow = FALSE)
lon_range <- lon_mat[1, ]
lat_range <- lat_mat[, 1]

# Cambiar LC_TIME para mostrar meses en inglés
old_locale <- Sys.getlocale("LC_TIME")
Sys.setlocale("LC_TIME", "C")

#===============================================================================
# 4. LOOP SOBRE CADA BIOMA
#===============================================================================
# Arreglos globales para unir resultados de cada bioma
status_matrix2 = status_matrix2_tot[,,193:264]
# global_BA_FireHarmonized_common_loyo <- array(NA, dim = dim(status_matrix2))
global_BA_FireHarmonized_common <- array(NA, dim = c(dim(status_matrix2_tot)[1],
dim(status_matrix2_tot)[2],length(1:72)))
global_BA_FireHarmonized_full    <- array(NA, dim = dim(status_matrix2_tot))

# status_matrix2 = status_matrix2_tot[,,193:264]
dim(global_BA_FireHarmonized_common)
dim(global_BA_FireHarmonized_full)
dim(status_matrix2)
# table_selected_all <- data.frame(
#   Bioma        = character(),
#   Month        = integer(),
#   Formula      = character(),
#   RFE_selected = character(),
#   stringsAsFactors = FALSE
# )
# Lista global para almacenar SHAP por bioma y mes
if (!exists("df_shap_global")) df_shap_global <- list()

biomas_unique <- unique(biomas_shp$cont_bm)

for (bioma in biomas_unique) {
  # bioma = "Africa-Mediterranean Forests, Woodlands & Scrub"  
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
  
  # BA_FireS3_crop_tot  <- BA_FireS3_tot[lon_idx,lat_idx,]
  # BA_Fire51_crop_tot  <- BA_Fire51_tot[lon_idx,lat_idx,]
  # 
  # BA_Fire51_crop_full   <- BA_Fire51_tot[lon_idx, lat_idx, ]
  # count_ActiveFire_crop_full <- count_ActiveFire_tot[lon_idx, lat_idx, ]
  
  # lat_mat_crop <- lat_mat[lat_idx, lon_idx]
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

  # st=status_crop[,,ind_common]
  # dim_crop <- dim(st)
  
  # Construir arrays para el período completo (usar recortes de totales)
  ntime_full <- dim(status_crop)[3]
  # obs_BA_Fire51_biome_full <- array(NA, dim = dim(status_crop))
  # obs_count_ActiveFire_biome_full <- array(NA, dim = dim(status_crop))
  
  
  dim_crop_tot <- dim(status_crop)
  # 
  # obs_BA_FireS3_biome_full <- array(NA, dim = status_crop)
  # obs_BA_Fire51_biome_full <- array(NA, dim = status_crop)
  # ---- Aplicar máscara a toda la serie ----
  n_time <- dim(status_crop)[3]
  status_masked <- array(NA, dim = c(length(lon_vec), length(lat_vec), n_time))
  
  for (t in 1:n_time) {
    slice <- status_crop[,,t]
    slice_masked <- ifelse(t(mask_matrix_sorted), slice, NA)
    status_masked[,,t] <- slice_masked
  }
  # dev.off()
  # ---- Visualización ejemplo de un mes (opcional) ----
  
  mask_final <- apply(status_masked, c(1,2), mean, na.rm = TRUE)

  mask_final_clean <- mask_final==1
  mask_final_clean[is.na(mask_final_clean)]=FALSE
  
  mask_final_clean <- apply(status_masked, c(1,2), function(x) any(x==1, na.rm = TRUE))
  # mask_final_clean <- mask_final==1
  mask_final_clean[is.na(mask_final_clean)]=FALSE
  # image(status_masked[,,2])
  # Guardar resultados en archivos RData (ajusta nombres de variables según convenga)
  # Sección de guardado (al final del loop por bioma)
  # load(paste0(output_dir_RData, "BA_", Modelo, "_", safe_biome_name, "_FireCCIS311.RData"))
  # load(paste0(output_dir_RData, "BA_", Modelo, "_", safe_biome_name, "_FireCCI51.RData"))

  # load(paste0(output_dir_RData, "BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Common_Loyo.RData"))
  load(paste0(output_dir_RData, "BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Common.RData"))
  load(paste0(output_dir_RData, "BA_", Modelo, "_", safe_biome_name, "_FireHarmonized_Full.RData"))
  
  # ==============================================================================
  # 7. UNIR RESULTADOS DEL BIOMA AL GLOBAL
  # ==============================================================================
  # mask_final[is.na(mask_final)] <- FALSE
  # dim(status_matrix2)
  # # Para el período común loyo
  # for(t in 1:dim(status_matrix2)[3]) {
  #   biome_values <- BA_FireHarmonized_common_loyo[,,t]    # Capa armonizada del bioma (común)
  #   tmp_global <- global_BA_FireHarmonized_common_loyo[lon_idx, lat_idx, t]
  #   tmp_global[ mask_final_clean ] <- biome_values[ mask_final_clean ]
  #   global_BA_FireHarmonized_common_loyo[lon_idx, lat_idx, t] <- tmp_global
  # }
  # # summary(as.vector(BA_FireHarmonized_common))
  # # dim(BA_FireHarmonized_common)
  # # Para el período común completo
  for(t in 1:dim(status_matrix2)[3]) {
    biome_values <- BA_FireHarmonized_common[,,t]    # Capa armonizada del bioma (común)
    tmp_global <- global_BA_FireHarmonized_common[lon_idx, lat_idx, t]
    tmp_global[ mask_final_clean ] <- biome_values[ mask_final_clean ]
    global_BA_FireHarmonized_common[lon_idx, lat_idx, t] <- tmp_global
  }
  # # Para el período completo
  for(t in 1:dim(status_matrix2_tot)[3]) {
    biome_values <- BA_FireHarmonized_full[,,t]    # Capa armonizada del bioma (completo)
    tmp_global <- global_BA_FireHarmonized_full[lon_idx, lat_idx, t]
    tmp_global[ mask_final_clean ] <- biome_values[ mask_final_clean ]
    global_BA_FireHarmonized_full[lon_idx, lat_idx, t] <- tmp_global
  }
  
  print(summary(as.vector(global_BA_FireHarmonized_common)))
  print(summary(as.vector(global_BA_FireHarmonized_common)))
  # print(summary(as.vector(global_BA_FireHarmonized_full)))
  cat("Resumen de BA_FireHarmonized (período común) para", safe_biome_name, ":\n")
  # print(summary(as.vector(BA_FireHarmonized_common)))
  capa_plot <- 8
  z_plot <- global_BA_FireHarmonized_full[,,capa_plot]
  
  if (sum(is.finite(z_plot)) > 0) {
    
    image.plot(
      lon,
      lat,
      z_plot,
      main = paste0("Global BA FireHarmonized - capa ", capa_plot)
    )
    
  } else {
    
    cat(
      "No se plotea global_BA_FireHarmonized_full[,,", capa_plot,
      "] porque no tiene valores finitos en este punto del loop.\n",
      sep = ""
    )
  }
  
}  # Fin loop por bioma
# === UNA VEZ FINALIZADO TODO EL BUCLE DE BIOMAS ===
image.plot(lon, lat, global_BA_FireHarmonized_full[,,2])
dev.off()
image.plot(lon, lat, global_BA_FireHarmonized_full[,,2])
dim(global_BA_FireHarmonized_full)
image.plot(lon, lat, global_BA_FireHarmonized_full[,,1])
summary(as.vector(global_BA_FireHarmonized_common))


dim(global_BA_FireHarmonized_common)
dim(global_BA_FireHarmonized_full)

# global_BA_FireHarmonized_common_loyo <- array(NA, dim = dim(status_matrix2))
# pp <- array(NA, dim = c(dim(global_BA_FireHarmonized_full)[1],
#                         dim(global_BA_FireHarmonized_full)[2],
#                         dim(global_BA_FireHarmonized_full)[3]))
# 

# ============================================================================
# 10. GUARDAR OBJETOS GLOBALES
# ============================================================================
# save(global_BA_FireHarmonized_common_loyo, lon_range, lat_range,
#      file = paste0(output_dir_RData, "BA_", Modelo, "global_BA_FireHarmonized_Common_Loyo.RData"))
save(global_BA_FireHarmonized_common, lon_range, lat_range,
file = paste0(output_dir_RData, "BA_", Modelo, "global_BA_FireHarmonized_Common.RData"))
save(global_BA_FireHarmonized_full, lon_range, lat_range, 
     file = paste0(output_dir_RData, "BA_", Modelo, "global_BA_FireHarmonized_Full.RData"))

#===============================================================================
# #perido común con loyo (global)
# ntime=72
# # 7. VISUALIZACIÓN DE LA SERIE DE TIEMPO
# dates <- seq(as.Date("2019-01-01"), by = "month", length.out = ntime)
# BA_Fire51<- BA_Fire51_tot[,,ind_common]
# 
# ba_total_Fire51 <- sapply(1:ntime, function(m) sum(BA_Fire51[,,m], na.rm = TRUE))
# ba_total_FireS3 <- sapply(1:ntime, function(m) sum(BA_FireS3[,,m], na.rm = TRUE))
# ba_total_harmonized_af <- sapply(1:ntime, function(m) sum(global_BA_FireHarmonized_common[,,m], na.rm = TRUE))
# 
# # Crear un data frame con las series temporales # Guardar el data frame
# series_temporales <- data.frame(tiempo = 1:ntime,
#                                 Ba_FireCCI51 = ba_total_Fire51,
#                                 Ba_FireCIIS311 = ba_total_FireS3,
#                                 Ba_Harmonized   = ba_total_harmonized_af)
# archivo_series <- paste0(output_dir_csv, "series_temporales_",Modelo,"_loyo_XGlobal.csv")
# write.csv(series_temporales, file = archivo_series, row.names = FALSE)
# 
# # Estadísticos para Harmonized vs S3
# y_true <- ba_total_FireS3
# y_pred_h <- ba_total_harmonized_af
# NME_h <- sum(abs(y_true - y_pred_h)) / sum(abs(y_true - mean(y_true)))
# bias_h   <- mean(y_pred_h - y_true)
# rmse_h   <- sqrt(mean((y_pred_h - y_true)^2))
# r2_h     <- summary(lm(y_pred_h ~ y_true))$r.squared
# cor_test_h <- cor.test(y_true, y_pred_h, method = "spearman")
# cor_h  <- as.numeric(cor_test_h$estimate)
# p_val_h    <- cor_test_h$p.value
# cor_signif_h <- ifelse(p_val_h < 0.05, "*", "ns")
# 
# # Estadísticos para F5 vs S3
# y_pred_f <- ba_total_Fire51
# bias_f   <- mean(y_pred_f - y_true)
# NME_f <- sum(abs(y_true - y_pred_f)) / sum(abs(y_true - mean(y_true)))
# rmse_f   <- sqrt(mean((y_pred_f - y_true)^2))
# r2_f     <- summary(lm(y_pred_f ~ y_true))$r.squared
# cor_test_f <- cor.test(y_true, y_pred_f, method = "spearman")
# cor_f  <- as.numeric(cor_test_f$estimate)
# p_val_f    <- cor_test_f$p.value
# cor_signif_f <- ifelse(p_val_f < 0.05, "*", "ns")
# 
# # Crear un data frame con los estadísticos y guardar
# estadisticas <- data.frame(
#   Modelo = c("F5_vs_S3", "Harmonized_vs_S3"), bias = c(bias_f, bias_h),
#   NME = c(NME_f, NME_h), RMSE = c(rmse_f, rmse_h),
#   R2 = c(r2_f, r2_h), Correlacion = c(cor_f, cor_h),
#   p_value = c(p_val_f, p_val_h),
#   Significancia = c(cor_signif_f, cor_signif_h))
# archivo_csv <- paste0(output_dir_csv, "Statistics_temporal_series_",Modelo,"_loyo_XGlobal.csv")
# write.csv(estadisticas, file = archivo_csv, row.names = FALSE)
# 
# dev.off()
# # # Crear y guardar la imagen
# jpeg_filename <- file.path(output_dir_plot, paste0("XGLOBAL_2019_2024_time_series_loyo.jpeg"))
# jpeg(filename = jpeg_filename, width = 1600, height = 800, res = 150)
# # Calcular el máximo del eje y
# max_y <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.5
# # Traza el plot base
# plot(dates, ba_total_FireS3, type = "l", lwd = 3, col = "blue",
#      xlab = "", ylab = "Burned Area (km²)", ylim = c(0, max_y), xaxt = "n")
# 
# # Definir los años únicos que abarca el vector de fechas
# years <- unique(format(dates, "%Y"))
# # Sombreado de verano (junio, julio y agosto) en gris claro
# for (year in years) {
#   summer_start <- as.Date(paste(year, "06", "01", sep = "-"))
#   summer_end   <- as.Date(paste(year, "08", "31", sep = "-"))
#   if (summer_end >= min(dates) && summer_start <= max(dates)) {
#     rect(max(min(dates), summer_start), 0,
#          min(max(dates), summer_end), max_y,
#          col = rgb(0.85, 0.85, 0.85, 0.5), border = NA)
#   }
# }
# 
# # Sombreado de invierno: enero y febrero, y diciembre, para cada año
# for (year in years) {
#   # Sombreado para enero y febrero del año actual
#   jan_start <- as.Date(paste(year, "01", "01", sep = "-"))
#   feb_end   <- as.Date(paste(year, "02", "28", sep = "-"))
#   if (feb_end >= min(dates) && jan_start <= max(dates)) {
#     rect(max(min(dates), jan_start), 0,
#          min(max(dates), feb_end), max_y,
#          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
#   }
# 
#   # Sombreado para diciembre del año actual
#   dec_start <- as.Date(paste(year, "12", "01", sep = "-"))
#   dec_end   <- as.Date(paste(year, "12", "31", sep = "-"))
#   if (dec_end >= min(dates) && dec_start <= max(dates)) {
#     rect(max(min(dates), dec_start), 0,
#          min(max(dates), dec_end), max_y,
#          col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
#   }
# }
# 
# # Agregar grilla y las líneas de las series para que queden sobre las áreas sombreadas
# grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
# lines(dates, ba_total_FireS3, col = "blue", lwd = 3)
# lines(dates, ba_total_Fire51, col = "orange", lwd = 3)
# lines(dates, ba_total_harmonized_af, col = "darkred", lwd = 3)
# axis.Date(1, at = dates, format = "%b", las = 2)
# 
# legend("top", legend = c("FireCCIS311", "FireCCI51", "Harmonized"),
#        col = c("blue", "orange", "darkred"), lwd = 3, ncol = 3, bty = "n", xpd = TRUE)
# 
# # Agregar textos de performance
# text_x_left <- max(dates) - (max(dates) - min(dates)) * 0.99
# text_y_left <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.22
# performance_text_left <- sprintf("Harmonized vs FireCCIS311:\nBias = %.0f\nNME = %.2f\nRMSE = %.0f\nR² = %.2f\nCor = %.2f %s",
#                                  bias_h, NME_h, rmse_h, r2_h, cor_h, cor_signif_h)
# text(text_x_left, text_y_left, performance_text_left, adj = c(0, 1), cex = 0.9, col = "darkred")
# 
# text_x_right <- max(dates) - (max(dates) - min(dates)) * 0.01
# text_y_right <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.22
# performance_text_right <- sprintf("FireCCI51 vs FireCCIS311:\nBias = %.0f\nNME = %.2f\nRMSE = %.0f\nR² = %.2f\nCor = %.2f %s",
#                                   bias_f, NME_f, rmse_f, r2_f, cor_f, cor_signif_f)
# text(text_x_right, text_y_right, performance_text_right, adj = c(1, 1), cex = 0.9, col = "darkred")
# dev.off()
#===============================================================================

ntime=72
# 7. VISUALIZACIÓN DE LA SERIE DE TIEMPO
dates <- seq(as.Date("2019-01-01"), by = "month", length.out = ntime)
BA_Fire51<- BA_Fire51_tot[,,ind_common]

ba_total_Fire51 <- sapply(1:ntime, function(m) sum(BA_Fire51[,,m], na.rm = TRUE))
ba_total_FireS3 <- sapply(1:ntime, function(m) sum(BA_FireS3[,,m], na.rm = TRUE))
ba_total_harmonized_af <- sapply(1:ntime, function(m) sum(global_BA_FireHarmonized_common[,,m], na.rm = TRUE))

# Crear un data frame con las series temporales # Guardar el data frame 
series_temporales <- data.frame(tiempo = 1:ntime, 
                                Ba_FireCCI51 = ba_total_Fire51,
                                Ba_FireCIIS311 = ba_total_FireS3,
                                Ba_Harmonized   = ba_total_harmonized_af)
archivo_series <- paste0(output_dir_csv, "series_temporales_",Modelo,"_XGlobal.csv")
write.csv(series_temporales, file = archivo_series, row.names = FALSE)

# Estadísticos para Harmonized vs S3
y_true <- ba_total_FireS3
y_pred_h <- ba_total_harmonized_af
NME_h <- sum(abs(y_true - y_pred_h)) / sum(abs(y_true - mean(y_true)))
bias_h   <- mean(y_pred_h - y_true)
rmse_h   <- sqrt(mean((y_pred_h - y_true)^2))
r2_h     <- summary(lm(y_pred_h ~ y_true))$r.squared
cor_test_h <- cor.test(y_true, y_pred_h, method = "spearman")
cor_h  <- as.numeric(cor_test_h$estimate)
p_val_h    <- cor_test_h$p.value
cor_signif_h <- ifelse(p_val_h < 0.05, "*", "ns")

# Estadísticos para F5 vs S3
y_pred_f <- ba_total_Fire51
bias_f   <- mean(y_pred_f - y_true)
NME_f <- sum(abs(y_true - y_pred_f)) / sum(abs(y_true - mean(y_true)))
rmse_f   <- sqrt(mean((y_pred_f - y_true)^2))
r2_f     <- summary(lm(y_pred_f ~ y_true))$r.squared
cor_test_f <- cor.test(y_true, y_pred_f, method = "spearman")
cor_f  <- as.numeric(cor_test_f$estimate)
p_val_f    <- cor_test_f$p.value
cor_signif_f <- ifelse(p_val_f < 0.05, "*", "ns")

# Crear un data frame con los estadísticos y guardar
estadisticas <- data.frame(
  Modelo = c("F5_vs_S3", "Harmonized_vs_S3"), bias = c(bias_f, bias_h),
  NME = c(NME_f, NME_h), RMSE = c(rmse_f, rmse_h),
  R2 = c(r2_f, r2_h), Correlacion = c(cor_f, cor_h),
  p_value = c(p_val_f, p_val_h), 
  Significancia = c(cor_signif_f, cor_signif_h))
archivo_csv <- paste0(output_dir_csv, "Statistics_temporal_series_",Modelo,"_XGlobal.csv")
write.csv(estadisticas, file = archivo_csv, row.names = FALSE)

dev.off()
# # Crear y guardar la imagen
jpeg_filename <- file.path(output_dir_plot, paste0("XGLOBAL_2019_2024_time_series.jpeg"))
jpeg(filename = jpeg_filename, width = 1600, height = 800, res = 150)
# Calcular el máximo del eje y
max_y <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.5
# Traza el plot base
plot(dates, ba_total_FireS3, type = "l", lwd = 3, col = "blue",
     xlab = "", ylab = "Burned Area (km²)", ylim = c(0, max_y), xaxt = "n")

# Definir los años únicos que abarca el vector de fechas
years <- unique(format(dates, "%Y"))
# Sombreado de verano (junio, julio y agosto) en gris claro
for (year in years) {
  summer_start <- as.Date(paste(year, "06", "01", sep = "-"))
  summer_end   <- as.Date(paste(year, "08", "31", sep = "-"))
  if (summer_end >= min(dates) && summer_start <= max(dates)) {
    rect(max(min(dates), summer_start), 0, 
         min(max(dates), summer_end), max_y, 
         col = rgb(0.85, 0.85, 0.85, 0.5), border = NA)
  }
}

# Sombreado de invierno: enero y febrero, y diciembre, para cada año
for (year in years) {
  # Sombreado para enero y febrero del año actual
  jan_start <- as.Date(paste(year, "01", "01", sep = "-"))
  feb_end   <- as.Date(paste(year, "02", "28", sep = "-"))
  if (feb_end >= min(dates) && jan_start <= max(dates)) {
    rect(max(min(dates), jan_start), 0, 
         min(max(dates), feb_end), max_y, 
         col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  }
  
  # Sombreado para diciembre del año actual
  dec_start <- as.Date(paste(year, "12", "01", sep = "-"))
  dec_end   <- as.Date(paste(year, "12", "31", sep = "-"))
  if (dec_end >= min(dates) && dec_start <= max(dates)) {
    rect(max(min(dates), dec_start), 0, 
         min(max(dates), dec_end), max_y, 
         col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)
  }
}

# Agregar grilla y las líneas de las series para que queden sobre las áreas sombreadas
grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
lines(dates, ba_total_FireS3, col = "blue", lwd = 3)
lines(dates, ba_total_Fire51, col = "orange", lwd = 3)
lines(dates, ba_total_harmonized_af, col = "darkred", lwd = 3)
axis.Date(1, at = dates, format = "%b", las = 2)

legend("top", legend = c("FireCCIS311", "FireCCI51", "Harmonized"),
       col = c("blue", "orange", "darkred"), lwd = 3, ncol = 3, bty = "n", xpd = TRUE)

# Agregar textos de performance
text_x_left <- max(dates) - (max(dates) - min(dates)) * 0.99
text_y_left <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.22
performance_text_left <- sprintf("Harmonized vs FireCCIS311:\nBias = %.0f\nNME = %.2f\nRMSE = %.0f\nR² = %.2f\nCor = %.2f %s", 
                                 bias_h, NME_h, rmse_h, r2_h, cor_h, cor_signif_h)
text(text_x_left, text_y_left, performance_text_left, adj = c(0, 1), cex = 0.9, col = "darkred")

text_x_right <- max(dates) - (max(dates) - min(dates)) * 0.01
text_y_right <- max(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af) * 1.22
performance_text_right <- sprintf("FireCCI51 vs FireCCIS311:\nBias = %.0f\nNME = %.2f\nRMSE = %.0f\nR² = %.2f\nCor = %.2f %s", 
                                  bias_f, NME_f, rmse_f, r2_f, cor_f, cor_signif_f)
text(text_x_right, text_y_right, performance_text_right, adj = c(1, 1), cex = 0.9, col = "darkred")
dev.off()


################################################################################

# Vector de fechas para el periodo total (2001-2024)
dates_tot <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
ntime_tot <- length(dates_tot)

BA_Fire51_tot[BA_Fire51_tot==0]=NA
BA_FireS3_tot[BA_FireS3_tot==0]=NA
global_BA_FireHarmonized_full[global_BA_FireHarmonized_full==0]=NA

# Series totales calculadas (suma de cada capa, con na.rm = TRUE)
ba_total_Fire51 <- sapply(1:ntime_tot, function(m) sum(BA_Fire51_tot[,,m], na.rm = TRUE))
ba_total_FireS3 <- sapply(1:ntime_tot, function(m) sum(BA_FireS3_tot[,,m], na.rm = TRUE))
ba_total_harmonized_af <- sapply(1:ntime_tot, function(m) sum(global_BA_FireHarmonized_full[,,m], na.rm = TRUE))

ba_total_FireS3[ba_total_FireS3==0]=NA

# Para FireCCIS311 (Reference) queremos usar solo datos de 2019-2024; 
# asignamos NA a los meses anteriores.
pre2019_idx <- which(dates_tot < as.Date("2019-01-01"))
# ba_total_FireS3_tot[pre2019_idx] <- NA
ba_total_FireS3[pre2019_idx] <- NA
# --- Definir periodos ---
# Periodo 1: 2001-2018 (para Harmonized)
period1_idx <- which(dates_tot < as.Date("2019-01-01"))
# Periodo 2: 2019-2024 (para Harmonized, FireCCI51 y Reference)
period2_idx <- which(dates_tot >= as.Date("2019-01-01"))

# --- Calcular medias ---
# Para Harmonized, se calculará en dos periodos:
mean_Harmonized_p1 <- mean(ba_total_harmonized_af[period1_idx], na.rm = TRUE)
mean_Harmonized_p2 <- mean(ba_total_harmonized_af[period2_idx], na.rm = TRUE)

# Para FireCCI51, se calcula la media solo para el periodo 2019-2024
mean_Fire51_p1 <- mean(ba_total_Fire51[period1_idx], na.rm = TRUE)
mean_Fire51_p2 <- mean(ba_total_Fire51[period2_idx], na.rm = TRUE)
# Para Reference (FireCCIS311), se calcula la media (ya que los NA se asignaron para 2001-2018)
mean_Reference <- mean(ba_total_FireS3, na.rm = TRUE)

series_temporales <- data.frame(
  tiempo = 1:ntime_tot,
  Ba_FireCCI51 = ba_total_Fire51,
  Ba_FireCIIS311 = ba_total_FireS3,
  Ba_Harmonized = ba_total_harmonized_af
)

archivo_series <- paste0(output_dir_csv, "series_temporales_", Modelo, "_XGlobal_tot.csv")
write.csv(series_temporales, file = archivo_series, row.names = FALSE)

dev.off()
# Crear y guardar el plot global total
jpeg_filename <- file.path(output_dir_plot, "XGLOBAL_2001_2024_time_series_tot.jpeg")
jpeg(filename = jpeg_filename, width = 1600, height = 800, res = 150)

# Calcular el máximo ignorando NA
max_y <- max(c(ba_total_FireS3, ba_total_Fire51, ba_total_harmonized_af), na.rm = TRUE) * 1.5

# Trazar el plot base usando dates_tot en todas las funciones
plot(dates_tot, ba_total_FireS3, type = "l", lwd = 1, col = "blue",
     xlab = "", ylab = "Burned Area (km²)",
     ylim = c(0, max_y), xaxt = "n")

lines(dates_tot, ba_total_Fire51, col = "orange", lwd = 1)
lines(dates_tot, ba_total_harmonized_af, col = "darkred", lwd = 1)



# --- Dibujar líneas horizontales que representen las medias ---
# Para FireCCI51 y Harmonized, la media se muestra para el periodo completo (2001-2024)
segments(x0 = as.Date("2001-01-01"), y0 = mean_Fire51_p1, 
         x1 = as.Date("2018-12-01"), y1 = mean_Fire51_p1,
         col = "orange", lwd = 2, lty = 2)
segments(x0 = as.Date("2019-01-01"), y0 = mean_Fire51_p2, 
         x1 = as.Date("2024-12-01"), y1 = mean_Fire51_p2,
         col = "orange", lwd = 2, lty = 2)


segments(x0 = as.Date("2001-01-01"), y0 = mean_Harmonized_p1, 
         x1 = as.Date("2018-12-31"), y1 = mean_Harmonized_p1,
         col = "darkred", lwd = 2, lty = 2)
segments(x0 = as.Date("2019-01-01"), y0 = mean_Harmonized_p2, 
         x1 = as.Date("2024-12-01"), y1 = mean_Harmonized_p2,
         col = "darkred", lwd = 2, lty = 3)

# Para FireCCIS311 (Reference), la línea se dibuja solo para 2019-2024:
segments(x0 = as.Date("2019-01-01"), y0 = mean_Reference, 
         x1 = as.Date("2024-12-01"), y1 = mean_Reference,
         col = "blue", lwd = 2, lty = 2)

# Configurar el eje X y la grilla:
axis.Date(1, at = seq(min(dates_tot), max(dates_tot), by = "3 months"),
          format = "%b-%Y", las = 2, cex.axis = 0.5)
grid(nx = NA, ny = NULL, col = "lightgray", lty = "dotted")

# Leyenda:
legend("topright", legend = c("FireCCIS311", 
                              "FireCCI51", 
                              "Harmonized"),
       col = c("blue", "orange", "darkred"),
       lwd = c(2, 2, 2), lty = c(2, 2, 2),
       ncol = 1, bty = "n", cex = 0.5, xpd = TRUE,
       seg.len = 1, text.width = NULL, inset = c(0, 0))

# Agregar texto con información de las medias en cada periodo:
usr <- par("usr")
performance_text <- sprintf(
  "Performance Comparisons:\n
   Mean (2001-2018): FireCCI51 (%.0f) | Harmonized (%.0f)\n
   Mean (2019-2024): FireCCIS311 (%.0f) | FireCCI51 (%.0f) | Harmonized (%.0f)",
  mean_Fire51_p1, mean_Harmonized_p1,
  mean_Reference, mean_Fire51_p2, mean_Harmonized_p2)

text(x = usr[1] + 50, y = usr[4],
     labels = performance_text,
     adj = c(0, 1), cex = 0.5, col = "black")

# dev.off()


# dev.off()

Sys.setlocale("LC_TIME", old_locale)



# 
# 
# # FireCCIS311 (disponible solo para el período común)
# load(paste0(dir_oss, "out-verifications-2019-2024_025/FireCCIS311_2019_2024_0.25degree-download.RData"))
# BA_FireS3 <- FireS3 / 1e6
# BA_FireS3[BA_FireS3 == 0] <- NA
# rm(FireS3)
# gc()
# dim(global_BA_FireHarmonized_common)
# 
# harmoni=global_BA_FireHarmonized_common
# harmoni[harmoni >= 0] <- 1
# # harmoni[is.na(harmoni)] <- 0
# 
# s3=BA_FireS3
# s3[s3 > 0] <- 1
# s3[is.na(s3)] <- 0
# 
# image.plot(lon, lat, s3[,,10])
# image.plot(lon, lat, harmoni[,,10], add=T)
# # Condición: s3 es 1 y harmoni no es 1
# diff_mask <- s3 == 1 & harmoni != 1
# # ¿Hay al menos un píxel que cumple esa condición?
# any(diff_mask, na.rm = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
