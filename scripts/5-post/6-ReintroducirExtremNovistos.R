
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
# Modelo <- "HBA-RF-2003-2024-CL60-R30"
Modelo <- "B1-MRBA60-2003-2024"
# Directorios de trabajo y salida
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
output_dir <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
# output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_RData <- paste0(output_dir, "/RData/")
output_dir_plot_rle <- paste0(output_dir, "/plot_ExtreNoHarmo/")

dirs <- c(output_dir, output_dir_csv, output_dir_plot, output_dir_RData, output_dir_plot_rle)
for (d in dirs) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# ---------------------------
# Cargar lon/lat y generar mallas
# ---------------------------
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")

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
ind_common   <- which(dates_full %in% dates_common)  # 2019–2024
months_vec   <- month(dates_full)

# ---------------------------
# Cargar FireCCIS311 (solo periodo común 2019-2024) y normalizar unidades
# ---------------------------
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

nrows <- dim(BA_Fire51_tot)[1]
ncols <- dim(BA_Fire51_tot)[2]
time_len <- dim(BA_Fire51_tot)[3]

# Alinear S3 al cubo completo con NA fuera del periodo común
BA_FireS3_tot <- array(NA_real_, dim = c(nrows, ncols, time_len))
BA_FireS3_tot[,,ind_common] <- BA_FireS3

# ---------------------------
# Cargar BA armonizado (pre-filtrado)
# ---------------------------
# ---------------------------
# Cargar BA armonizado (pre-filtrado)
# ---------------------------
load(paste0(output_dir_RData, "/BA_B1-MRBA60-2003-2024global_BA_FireHarmonized_Full.RData"))
BA_harmonised <- global_BA_FireHarmonized_full
image.plot(lon,lat, BA_harmonised[,,1])
rm(global_BA_FireHarmonized_full); gc()
# load(paste0(output_dir_RData, "BA_M3_RF-PROB-2003-2024-CL60-R30global_BA_FireHarmonized_Full.RData"))
# BA_harmonised <- global_BA_FireHarmonized_full
BA_harmonised[BA_harmonised == 0] <- NA
# dim(BA_harmonised)
# BA_harmonised[is.na(BA_harmonised)] <- 0
image.plot(lon, lat, BA_harmonised[,,12])
dev.off()
# load(paste0(output_dir_RData, "BA_M3_RF-PROB-2003-2024-CL60-R30global_BA_FireHarmonized_Full.RData"))
# BA_harmonised <- global_BA_FireHarmonized_full
# BA_harmonised[BA_harmonised == 0] <- NA
# # BA_harmonised[is.na(BA_harmonised)] <- 0
# image.plot(lon, lat, BA_harmonised[,,12])
# rm(BA_final); gc()
# dev.off()

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
# load(paste0(dir_oss, "status_matrix2_tot.RData"))
# status_matrix2_tot=status_matrix2_tot[,,25:264]
BA_Fire51_tot[status_matrix2_tot == 1 & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[status_matrix2_tot == 1 & is.na(BA_FireS3_tot)] <- 0
BA_harmonised[status_matrix2_tot == 1 & is.na(BA_harmonised)] <- 0
gc()




# ============================================================================
# biomas_shp <- st_read("/mnt/disco6tb/FireCCI60/data_025/continental-biomes_dinerstein_V9/continental-biomes_dinerstein_V10.shp")
# biomas_shp <- st_transform(biomas_shp, crs = 4326)

biomas_shp <- st_read(file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"))
biomas_shp <- st_transform(biomas_shp, crs = 4326)

biomas_unique <- unique(biomas_shp$cont_bm)





# === Rutas para leer XLSX de topes y guardar resultados ===
suppressPackageStartupMessages({ library(readxl) })  # por si no estaba cargado

# Los XLSX están en la carpeta /csv del modelo activo:
dir_xlsx_main <- output_dir_csv          # p.ej. ".../FireCCI60-2003-2024-CL30-R30/csv/"
dir_xlsx_alt  <- dir_xlsx_main            # por si has adjuntado algún xlsx aquí

# Carpeta de salida para máscaras y agregados
out_dir_exceed <- file.path(output_dir_RData, "threshold_exceedance_F51")
dir.create(out_dir_exceed, recursive = TRUE, showWarnings = FALSE)

# (opcional) imprime para verificar
cat("Leyendo topes desde:", dir_xlsx_main, "\n")
cat("Guardando salidas en:", out_dir_exceed, "\n")

# Fechas y vectores mes/año (2003-01 a 2024-12)
date_seq <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
stopifnot(length(date_seq) == dim(BA_Fire51_tot)[3])
moy_vec  <- as.integer(format(date_seq, "%m")) # 1..12
yr_vec   <- as.integer(format(date_seq, "%Y"))

# (Opcional) contenedor en memoria
exceed_results <- list()
# --- Fechas (2003-01 a 2024-12) ---
date_seq <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
stopifnot(length(date_seq) == dim(BA_Fire51_tot)[3])
moy_vec  <- as.integer(format(date_seq, "%m"))   # 1..12
yr_vec   <- as.integer(format(date_seq, "%Y"))

# --- Acumuladores GLOBALS (mismo tamaño que BA_Fire51_tot) ---
BA_Fire51_aboveTope_global     <- array(NA_real_, dim = dim(BA_Fire51_tot))   # valores
Mask_Fire51_aboveTope_global   <- array(FALSE,     dim = dim(BA_Fire51_tot))  # binaria

cat("Leyendo topes desde:", dir_xlsx_main, "\n")
cat("Guardando salidas globales en:", output_dir_RData, "\n")

# =========================
# LOOP POR BIOMAS (tu lógica)
# =========================
for (bioma in biomas_unique) {
  cat("\nProcesando bioma:", bioma, "\n")
  safe_biome_name <- gsub(" ", "_", bioma)
  safe_biome_name <- gsub("[^[:alnum:]_]", "", safe_biome_name)
  
  # --- Selección bioma y bbox (manteniendo tu método) ---
  bioma_sel <- biomas_shp %>% dplyr::filter(cont_bm == bioma)
  bbox <- st_bbox(bioma_sel)
  
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  if (length(lon_idx) == 0 || length(lat_idx) == 0) { cat("Sin celdas para", bioma, "\n"); next }
  
  lon_mat_crop <- lon_mat[lat_idx, lon_idx]
  lat_mat_crop <- lat_mat[lat_idx, lon_idx]
  lon_vec <- sort(unique(as.vector(lon_mat_crop)))
  lat_vec <- sort(unique(as.vector(lat_mat_crop)))
  
  status_crop <- status_matrix2_tot[lon_idx, lat_idx, ]
  
  # ---- Crear objeto espacial y asignar bioma (tu flujo) ----
  grid_points_crop <- st_as_sf(data.frame(
    lon = as.vector(lon_mat_crop),
    lat = as.vector(lat_mat_crop)),
    coords = c("lon", "lat"), crs = 4326
  )
  inter <- st_intersects(grid_points_crop, biomas_shp)
  bioma_asignado <- sapply(inter, function(i) {
    if (length(i) == 0) return("Ninguno")
    return(biomas_shp$cont_bm[i[1]])
  })
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  presencia_en_serie <- apply(status_crop, c(1,2), function(x) any(x == 1, na.rm = TRUE))
  grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie))
  
  puntos_sin_bioma_con_presencia <- grid_points_biomas %>% dplyr::filter(bioma_final == "Ninguno" & presencia_serie)
  puntos_con_bioma <- grid_points_biomas %>% dplyr::filter(bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    nearest_idx <- st_nearest_feature(puntos_sin_bioma_con_presencia, puntos_con_bioma)
    grid_points_biomas$bioma_final[
      which(grid_points_biomas$bioma_final == "Ninguno" & grid_points_biomas$presencia_serie)
    ] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  # ---- Máscara final del bioma (tu método) ----
  mask <- grid_points_biomas$bioma_final == bioma
  mask_matrix <- matrix(mask, nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE)
  mask_matrix_sorted <- mask_matrix[order(lat_vec), order(lon_vec)]
  
  n_time <- dim(status_crop)[3]
  status_masked <- array(NA, dim = c(length(lon_vec), length(lat_vec), n_time))
  for (t in 1:n_time) {
    slice <- status_crop[,,t]
    slice_masked <- ifelse(t(mask_matrix_sorted), slice, NA)
    status_masked[,,t] <- slice_masked
  }
  mask_final <- apply(status_masked, c(1,2), mean, na.rm = TRUE)
  mask_final_clean <- mask_final == 1
  mask_final_clean[is.na(mask_final_clean)] <- FALSE
  mask_final_sum <- (t(mask_matrix) == TRUE | mask_final_clean == TRUE)  # [lon, lat] lógico
  
  # --- Recorte BA de F51 (mismo bbox) ---
  BA_Fire51_crop_tot <- BA_Fire51_tot[lon_idx, lat_idx, ]
 
  
  file_match <- NA_character_
  
  find_tope_file <- function(dir_topes, bioma) {
    
    if (!dir.exists(dir_topes)) return(NA_character_)
    
    all_files <- list.files(
      dir_topes,
      pattern = "^Max_Tope_COMMON_ByMonth_.*\\.(csv|xlsx)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    
    if (length(all_files) == 0) return(NA_character_)
    
    biome_clean <- tolower(gsub("[^[:alnum:]]+", "", bioma))
    
    for (ff in all_files) {
      
      base <- basename(ff)
      base <- sub("^Max_Tope_COMMON_ByMonth_", "", base)
      base <- sub("\\.(csv|xlsx)$", "", base, ignore.case = TRUE)
      base_clean <- tolower(gsub("[^[:alnum:]]+", "", base))
      
      if (
        grepl(biome_clean, base_clean, fixed = TRUE) ||
        grepl(base_clean, biome_clean, fixed = TRUE)
      ) {
        return(ff)
      }
    }
    
    return(NA_character_)
  }
  
  # Buscar en carpeta principal
  file_match <- find_tope_file(dir_xlsx_main, bioma)

  
  # Si no se encuentra, informar y continuar con el siguiente bioma
  if (is.na(file_match)) {
    
    cat("No se encontró archivo de topes CSV/XLSX para:", bioma, " — se omite.\n")
    
    cat("\nArchivos disponibles en dir_xlsx_main:\n")
    print(list.files(
      dir_xlsx_main,
      pattern = "^Max_Tope_COMMON_ByMonth_.*",
      full.names = FALSE
    ))
    
    if (!identical(dir_xlsx_main, dir_xlsx_alt) && dir.exists(dir_xlsx_alt)) {
      cat("\nArchivos disponibles en dir_xlsx_alt:\n")
      print(list.files(
        dir_xlsx_alt,
        pattern = "^Max_Tope_COMMON_ByMonth_.*",
        full.names = FALSE
      ))
    }
    
    next
  }
  
  cat("Topes:", basename(file_match), "\n")
  
  # =========================
  # LEER TABLA SEGÚN EXTENSIÓN
  # =========================
  
  ext_file <- tolower(tools::file_ext(file_match))
  
  if (ext_file == "xlsx") {
    
    suppressPackageStartupMessages({
      library(readxl)
    })
    
    tb <- readxl::read_xlsx(file_match)
    
  } else if (ext_file == "csv") {
    
    tb <- read.csv(
      file_match,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else {
    
    stop("Extensión no soportada para archivo de topes: ", file_match)
  }
  
  # =========================
  # DETECTAR COLUMNA DE TOPE Y MES
  # =========================
  
  names_clean <- tolower(trimws(names(tb)))
  
  # Usar el máximo de FireCCI51 en el periodo común
  tope_col_idx <- which(names_clean == "max_f51_in_common")
  
  # Si en algún momento quieres usar FireCCIS311/S3, cambia por:
  # tope_col_idx <- which(names_clean == "max_f3_in_common")
  
  # Si tienes una columna combinada, podrías usar:
  # tope_col_idx <- which(names_clean %in% c("tope_max_f51_f3", "max_tope_common"))
  
  if (length(tope_col_idx) == 0) {
    stop(
      "No se encontró la columna 'Max_F51_in_COMMON' en: ",
      file_match,
      "\nColumnas disponibles: ",
      paste(names(tb), collapse = ", ")
    )
  }
  
  tope_vec_raw <- suppressWarnings(as.numeric(tb[[tope_col_idx[1]]]))
  
  mois_idx <- which(names_clean %in% c("mes", "month", "moy", "month_of_year"))
  
  if (length(mois_idx) > 0) {
    
    meses_tb <- suppressWarnings(as.integer(tb[[mois_idx[1]]]))
    
  } else {
    
    if (length(tope_vec_raw) < 12) {
      stop("La tabla de topes no tiene al menos 12 filas mensuales: ", file_match)
    }
    
    meses_tb <- 1:12
    tope_vec_raw <- tope_vec_raw[1:12]
  }
  
  # =========================
  # VECTOR DE 12 TOPES MENSUALES
  # =========================
  
  tope_by_month <- rep(NA_real_, 12)
  
  for (mm in 1:12) {
    idx_mm <- which(meses_tb == mm)[1]
    
    if (!is.na(idx_mm)) {
      tope_by_month[mm] <- tope_vec_raw[idx_mm]
    }
  }
  
  if (any(is.na(tope_by_month))) {
    stop(
      "Faltan topes para algún mes en: ",
      file_match,
      "\nMeses encontrados: ",
      paste(meses_tb, collapse = ", ")
    )
  }
  
  cat("Topes mensuales cargados correctamente:\n")
  print(data.frame(
    month = 1:12,
    threshold = tope_by_month
  ))
  # =========================
  # APLICAR TOPES Y ACUMULAR EN LOS GLOBALES
  # =========================
  for (t in 1:dim(BA_Fire51_crop_tot)[3]) {
    mm <- moy_vec[t]
    thr <- tope_by_month[mm]
    
    slice_ba <- BA_Fire51_crop_tot[,,t]
    exceed <- (slice_ba > thr)
    exceed[is.na(exceed)] <- FALSE
    
    # limitar al bioma
    exceed <- exceed & mask_final_sum
    
    # actualizar máscara global (OR)
    Mask_Fire51_aboveTope_global[lon_idx, lat_idx, t] <- Mask_Fire51_aboveTope_global[lon_idx, lat_idx, t] | exceed
    
    # actualizar valores globales solo donde excede; mantener lo ya escrito en otras iteraciones
    prev_vals <- BA_Fire51_aboveTope_global[lon_idx, lat_idx, t]
    if (any(exceed)) prev_vals[exceed] <- slice_ba[exceed]
    BA_Fire51_aboveTope_global[lon_idx, lat_idx, t] <- prev_vals
  }
}

# image.plot(lon, lat, apply(BA_Fire51_aboveTope_global, c(1,2), mean, na.rm=T))
# summary(as.vector(BA_Fire51_aboveTope_global))
dev.off()
# =========================
# GUARDAR LOS DOS RDATA GLOBALES
# =========================
file_vals  <- file.path(output_dir_RData, paste0("F51_aboveTope_GLOBAL_", Modelo, "_2003_2024.RData"))
file_mask  <- file.path(output_dir_RData, paste0("F51_maskAboveTope_GLOBAL_", Modelo, "_2003_2024.RData"))

save(BA_Fire51_aboveTope_global, file = file_vals)
save(Mask_Fire51_aboveTope_global, file = file_mask)

cat("\nGuardado global de valores  :", file_vals,  "\n")
cat("Guardado global de máscara :", file_mask,  "\n")

dim(BA_Fire51_aboveTope_global)



load(paste0(output_dir_RData, "/BA_B1-MRBA60-2003-2024global_BA_FireHarmonized_Full.RData"))
BA_harmonised <- global_BA_FireHarmonized_full
dim(BA_harmonised)
image.plot(lon,lat, BA_harmonised[,,1])
rm(global_BA_FireHarmonized_full); gc()
# image(Mask_Fire51_aboveTope_global[,,1])
dev.off()
# Índices de tiempo a aplicar (2003-2018)
time_idx_apply <- 1:192

# Por si la máscara está como 0/1 en lugar de TRUE/FALSE
Mask_Fire51_aboveTope_global <- (Mask_Fire51_aboveTope_global == 1) | (Mask_Fire51_aboveTope_global == TRUE)

cat("Aplicando sustitución en 2003–2018 (t = 1:192)...\n")

for (t in time_idx_apply) {
  m <- Mask_Fire51_aboveTope_global[,,t]
  if (!any(m, na.rm = TRUE)) next
  # Sustituir solo donde la máscara es TRUE
  BA_harmonised_slice <- BA_harmonised[,,t]
  F51_vals_slice      <- BA_Fire51_aboveTope_global[,,t]
  BA_harmonised_slice[m] <- F51_vals_slice[m]
  BA_harmonised[,,t] <- BA_harmonised_slice
}

dim(BA_harmonised)
# Guardar resultado corregido
out_file_corr <- file.path(output_dir_RData, paste0("BA_harmonised_correctedByF51Tope_", Modelo, ".RData"))
save(BA_harmonised, file = out_file_corr)
cat("Guardado:\n  ", out_file_corr, "\n")



# dim(Mask_Fire51_aboveTope_global)
# image(Mask_Fire51_aboveTope_global[,,1])



# --- Porcentajes de restauración ---

# Índices 2003–2018 (ya los tienes)
time_idx_apply <- 1:192

# Conteos por tiempo
restored_counts <- sapply(time_idx_apply, function(t) {
  sum(Mask_Fire51_aboveTope_global[,,t], na.rm = TRUE)
})

# Denominador por tiempo: celdas válidas (no NA) en BA_harmonised ese mes
total_counts <- sapply(time_idx_apply, function(t) {
  sum(!is.na(BA_harmonised[,,t]))
})

# % restaurado por mes
pct_by_time <- 100 * restored_counts / total_counts

# % restaurado global espacio×tiempo (2003–2018)
overall_pct <- 100 * sum(restored_counts) / sum(total_counts)

cat(sprintf("%% restaurado global 2003–2018: %.3f%%\n", overall_pct))

# (Opcional) guardar serie mensual
pct_df <- data.frame(
  date = seq(as.Date("2003-01-01"), by = "month", length.out = length(time_idx_apply)),
  pct_restaurado = pct_by_time
)
write.csv(pct_df, file.path(output_dir_RData, paste0("pct_restaurado_mensual_2003_2018_", Modelo, ".csv")),
          row.names = FALSE)

# --- % de celdas del grid que ALGUNA VEZ fueron restauradas (2003–2018) ---

# Celdas con restauración alguna vez en 2003–2018
any_restored <- apply(Mask_Fire51_aboveTope_global[,,time_idx_apply], c(1,2), any)

# Celdas válidas al menos en un mes (evita contar mar/NA)
valid_any <- apply(!is.na(BA_harmonised[,,time_idx_apply]), c(1,2), any)

pct_cells_ever <- 100 * sum(any_restored & valid_any, na.rm = TRUE) / sum(valid_any, na.rm = TRUE)
cat(sprintf("%% de celdas válidas que alguna vez fueron restauradas (2003–2018): %.3f%%\n", pct_cells_ever))

# --- Frecuencia espacial de restauración (% de meses restaurada por celda) ---
# Promedio temporal de la máscara (TRUE=1/FALSE=0) por celda
freq_restored <- apply(Mask_Fire51_aboveTope_global[,,time_idx_apply], c(1,2), function(x) mean(x, na.rm = TRUE))
freq_restored_pct <- 100 * freq_restored
# freq_restored_pct
# (Opcional) visualización rápida
# image.plot(lon, lat, freq_restored_pct, main = "% de meses restaurada (2003–2018)")

# # Guardar mapas y métricas
# save(pct_by_time, overall_pct, any_restored, pct_cells_ever, freq_restored_pct,
#      file = file.path(output_dir_RData, paste0("resumen_restauracion_2003_2018_", Modelo, ".RData")))
dev.off()

