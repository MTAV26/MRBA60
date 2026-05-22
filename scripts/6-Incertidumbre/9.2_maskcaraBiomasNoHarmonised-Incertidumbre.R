

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

# ---------------------------
# Modelo <- "HBA-RF-2003-2022-CL60-R30"
# Modelo <- "FireCCI60-2003-2022-CL30-R30-3"
# # Directorios de trabajo y salida
# dir_oss <- '/mnt/disco6tb/FireCCI60/data_025/'
# output_dir <- paste0("/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/", Modelo)
Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
output_dir <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
output_dir_RData <- paste0(output_dir, "/RData/")
output_dir_plot_rle <- paste0(output_dir, "/plot_QM001M/")


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
ind_common   <- which(dates_full %in% dates_common)  # 2019–2022
months_vec   <- month(dates_full)

# ---------------------------
# Cargar FireCCIS311 (solo periodo común 2019-2022) y normalizar unidades
# ---------------------------

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


nrows <- dim(BA_Fire51_tot)[1]
ncols <- dim(BA_Fire51_tot)[2]
time_len <- dim(BA_Fire51_tot)[3]

# Alinear S3 al cubo completo con NA fuera del periodo común
BA_FireS3_tot <- array(NA_real_, dim = c(nrows, ncols, time_len))
BA_FireS3_tot[,,ind_common] <- BA_FireS3

# ---------------------------
# Cargar BA armonizado (pre-filtrado)
# ---------------------------
# Guardar resultado corregido

load(paste0(output_dir_RData, "/BA_harmonised_correctedByF51Tope_B1-MRBA60-2003-2024.RData"))
BA_harmonised[BA_harmonised == 0] <- NA
# BA_harmonised[is.na(BA_harmonised)] <- 0
image.plot(lon, lat, BA_harmonised[,,12])
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
gc()

library(dplyr)
library(ggplot2)


# =========================
# TIEMPO / CALENDARIO
# =========================
# Asumimos mensual desde 2003-01-01 (ajusta si tu NetCDF trae otras fechas)
start_date <- as.Date("2003-01-01")
nt <- dim(BA_Fire51_tot)[3]
dates_all <- seq.Date(from = start_date, by = "month", length.out = nt)
mons_all  <- month(dates_all)



# ============================================================================

biomas_shp <- st_read(file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"))
biomas_shp <- st_transform(biomas_shp, crs = 4326)
biomas_unique <- unique(biomas_shp$cont_bm)
# ============================================================
# MÁSCARA GLOBAL (TRUE/FALSE) POR BIOMAS Y MESES — ESTILO TU SCRIPT
# ============================================================

# --- meses por bioma (1=Ene ... 12=Dic) ---
selected_biome_months <- list(
  "Australia-Tropical Broadleaf Forests"                    = c(2, 3, 4),               # Feb, Mar, Abr
  "Eurasia-Tundra"                                          = c(1, 2, 12),             # Ene, Nov, Dic
  "Europe-Boreal Forests/Taiga"                             = c(1, 2, 12),          # Ene, Feb, Nov, Dic
  "North America-Boreal Forests/Taiga"                      = c(2),                     # Feb
  "North America-Mediterranean Forests, Woodlands & Scrub"  = c(1),                     # Ene
  "North America-Tundra"                                    = c(1, 2, 3, 4, 10, 11, 12),# Ene, Feb, Mar, Abr, Oct, Nov, Dic
  "South America-Temperate Broadleaf & Mixed Forests"       = c(7,8,9)                      # Mar
)

# --- eje temporal a usar (meses) ---
if (exists("months_vec")) {
  months_global <- months_vec
} else if (exists("mons_all")) {
  months_global <- mons_all
} else {
  stop("No encuentro months_vec ni mons_all. Define uno de los dos con los meses (1..12) del eje temporal.")
}

# --- máscara global inicial ---
mask_truefalse <- array(FALSE, dim = dim(BA_Fire51_tot))  # [lon, lat, time]
# 
# # ============================================================================
# biomas_shp <- st_read("/mnt/disco6tb/FireCCI60/data_025/continental-biomes_dinerstein_V9/continental-biomes_dinerstein_V10.shp")
# biomas_shp <- st_transform(biomas_shp, crs = 4326)



# --- biomas disponibles ---
# biomas_unique <- sort(unique(biomas_shp$cont_bm))

# =========================
# LOOP POR BIOMAS
# =========================
for (bioma in biomas_unique) {
  # si el bioma no está en la lista de meses seleccionados, lo saltamos (será FALSE)
  if (!(bioma %in% names(selected_biome_months))) next
  
  cat("\nProcesando bioma:", bioma, "\n")
  safe_biome_name <- gsub(" ", "_", bioma)
  safe_biome_name <- gsub("[^[:alnum:]_]", "", safe_biome_name)
  
  # --- Selección bioma y bbox ---
  bioma_sel <- biomas_shp %>% dplyr::filter(cont_bm == bioma)
  if (nrow(bioma_sel) == 0) { cat("  * No se encontró en biomas_shp\n"); next }
  bbox <- st_bbox(bioma_sel)
  
  # --- Índices dentro del bbox ---
  lon_idx <- which(lon >= bbox["xmin"] & lon <= bbox["xmax"])
  lat_idx <- which(lat >= bbox["ymin"] & lat <= bbox["ymax"])
  if (length(lon_idx) == 0 || length(lat_idx) == 0) { cat("  * BBox sin celdas\n"); next }
  
  # --- Recortes de coordenadas ---
  lon_mat_crop <- lon_mat[lat_idx, lon_idx, drop = FALSE]  # [lat,lon]
  lat_mat_crop <- lat_mat[lat_idx, lon_idx, drop = FALSE]  # [lat,lon]
  lon_vec_crop <- as.vector(lon_mat_crop)
  lat_vec_crop <- as.vector(lat_mat_crop)
  lon_vec <- sort(unique(lon_vec_crop))
  lat_vec <- sort(unique(lat_vec_crop))
  
  # --- status (tu lógica de “celdas con info” ya creada antes) ---
  status_crop <- status_matrix2_tot[lon_idx, lat_idx, , drop = FALSE]  # [lon,lat,time]
  
  # ---- Crear objeto espacial de puntos del recorte ----
  grid_points_crop <- st_as_sf(
    data.frame(lon = lon_vec_crop, lat = lat_vec_crop),
    coords = c("lon", "lat"), crs = 4326
  )
  
  # ---- Asignación de biomas (primer match) ----
  inter <- st_intersects(grid_points_crop, biomas_shp)
  bioma_asignado <- sapply(inter, function(i) {
    if (length(i) == 0) return("Ninguno")
    return(biomas_shp$cont_bm[i[1]])
  })
  grid_points_biomas <- grid_points_crop
  grid_points_biomas$bioma_final <- bioma_asignado
  
  # ---- Reasignar puntos “Ninguno” con presencia en la serie ----
  # presencia con tu status (1 = hay información en algún producto/arm)
  presencia_en_serie <- apply(status_crop, c(1,2), function(x) any(x == 1, na.rm = TRUE))  # [lon,lat]
  grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie)) # a orden [lat,lon] vectorizado por columnas
  
  puntos_sin_bioma_con_presencia <- grid_points_biomas %>%
    dplyr::filter(bioma_final == "Ninguno" & presencia_serie)
  
  puntos_con_bioma <- grid_points_biomas %>%
    dplyr::filter(bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    nearest_idx <- st_nearest_feature(puntos_sin_bioma_con_presencia, puntos_con_bioma)
    grid_points_biomas$bioma_final[
      which(grid_points_biomas$bioma_final == "Ninguno" & grid_points_biomas$presencia_serie)
    ] <- puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  # ---- Máscara espacial del bioma seleccionado (matriz lógica) ----
  mask <- grid_points_biomas$bioma_final == bioma
  mask_matrix <- matrix(mask, nrow = length(lat_vec), ncol = length(lon_vec), byrow = FALSE) # [lat,lon] en orden “geom”
  mask_matrix_sorted <- mask_matrix[order(lat_vec), order(lon_vec)]
  mask_final_clean <- t(mask_matrix_sorted)  # [lon,lat] TRUE dentro del bioma
  
  # --- Meses a marcar para este bioma ---
  meses_ok <- selected_biome_months[[bioma]]
  t_idx <- which(months_global %in% meses_ok)
  if (length(t_idx) == 0) next
  
  # --- Volcar al cubo global (tu mismo patrón) ---
  for (tt in t_idx) {
    tmp_global <- mask_truefalse[lon_idx, lat_idx, tt, drop = FALSE][,,1]
    tmp_global[mask_final_clean] <- TRUE
    mask_truefalse[lon_idx, lat_idx, tt] <- tmp_global
  }
  
  image(mask_truefalse[,,2])
}

cat("\nMáscara global creada.\n",
    "Dimensiones:", paste(dim(mask_truefalse), collapse = " x "), "\n",
    "Proporción de TRUE:", round(mean(mask_truefalse), 6), "\n")

# (Opcional) Guardar en disco:
save(mask_truefalse, file = file.path(output_dir_RData, paste0("mask_truefalse_biomas_meses_", Modelo, ".RData")))

# (Opcional) Verificación rápida:
# k <- 1L  # índice temporal a mirar
# image.plot(lon, lat, t(mask_truefalse[,,k]), main = paste("Máscara TRUE/FALSE - t =", k))





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
  library(scales)
})

# installed.packages("qmap")
# ---------------------------
# Configuración y rutas
# ---------------------------
# Modelo <- "HBA-RF-2003-2022-CL60-R30"
Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
output_dir <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
output_dir_RData <- paste0(output_dir, "/RData/")
output_dir_plot_rle <- paste0(output_dir, "/plot_QM001M/")


# ---------------------------
# Cargar lon/lat y generar mallas
# ---------------------------
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")


load(file.path(output_dir_RData, "BA_Incertidumbre_FireCCI51.RData"))
incert_Fire51=global_rmse_abs_full
rm(global_rmse_abs_full)
gc()

load(file.path(output_dir_RData, "BA_Incertidumbre_HARMONISED_abs.RData"))
incert_Harmo=global_rmse_abs_full
rm(global_rmse_abs_full)
gc()

# file_mask  <- file.path(output_dir_RData, paste0("F51_maskAboveTope_GLOBAL_", Modelo, "_2003_2022.RData"))
load(paste0(output_dir_RData, "/F51_maskAboveTope_GLOBAL_", Modelo, "_2003_2024.RData"))
dim(Mask_Fire51_aboveTope_global)
image(Mask_Fire51_aboveTope_global[,,2])
Mask_Fire51_aboveTope_global <- (Mask_Fire51_aboveTope_global == 1) | (Mask_Fire51_aboveTope_global == TRUE)
Mask_Fire51_aboveTope_global[is.na(Mask_Fire51_aboveTope_global)] <- FALSE
Mask_Fire51_aboveTope_global[,,1]
dim(Mask_Fire51_aboveTope_global)

load(paste0(output_dir_RData,"mask_truefalse_biomas_meses_", Modelo, ".RData"))
image(mask_truefalse[,,1])
mask_truefalse[,,1]
dim(mask_truefalse)

dim(Mask_Fire51_aboveTope_global)
dim(incert_Harmo)
dim(incert_Fire51)
image(Mask_Fire51_aboveTope_global[,,1])

Mask_union_global<-Mask_Fire51_aboveTope_global | mask_truefalse
image(Mask_union_global[,,2])

incert_Fire60=incert_Harmo
nt <- dim(incert_Fire60)[3]

for (t in seq_len(nt)) {
  idx <- (Mask_union_global[,,t] == 1)
  h <- incert_Fire60[,,t]
  h[idx] <- h[idx] + incert_Fire51[,,t][idx]
  incert_Fire60[,,t] <- h
}
# 
# outfile <- file.path(output_dir_RData, "BA_Incertidumbre_FireCCI60.RData")
# save(incert_Fire60, file = outfile)
image.plot(incert_Fire60[,,1])



# ============================================================
# LECTURA DE CAPAS Y REEMPLAZO SEGÚN MÁSCARA
# ============================================================

# # --- 1) Leer máscara ---
# mask_file <- file.path(output_dir_RData, paste0("mask_truefalse_biomas_meses_", Modelo, ".RData"))
# load(mask_file)  # carga: mask_truefalse

# Asegurar tipo lógico
if (!is.logical(mask_truefalse)) {
  mask_truefalse <- Mask_union_global != 0
}
image.plot(Mask_union_global[,,4])
out_replaced <- file.path(output_dir_RData, "MASK_NOHARMONISED.RData")
save(Mask_union_global, file = out_replaced)

dim(mask_truefalse)

# # --- 3) Leer Incertidumbre FireCCI51 (como indicas) ---
# load(file.path(output_dir_RData, "BA_Incertidumbre_FireCCI51.RData"))
# if (!exists("global_rmse_abs_full")) stop("BA_Incertidumbre_FireCCI51.RData no contiene 'global_rmse_abs_full'.")
# incert_Fire51 <- global_rmse_abs_full
# rm(global_rmse_abs_full); gc()

# --- 4) Comprobaciones de dimensiones ---
d60  <- dim(incert_Fire60)
d51  <- dim(incert_Fire51)
dmsk <- dim(mask_truefalse)
# if (is.null(d60) || is.null(d51) || is.null(dmsk)) stop("Alguno de los objetos no es array/matriz con dim definidas.")
# 
# if (!all(d60 == d51))  stop("Dimensiones de incert_Fire60 y incert_Fire51 no coinciden: ",
#                             paste(d60, collapse="x"), " vs ", paste(d51, collapse="x"))
# if (!all(d60 == dmsk)) stop("Dimensiones de incert_Fire60 y mask_truefalse no coinciden: ",
#                             paste(d60, collapse="x"), " vs ", paste(dmsk, collapse="x"))

# --- 5) Reemplazo vectorizado donde mask == TRUE ---
n_total    <- length(mask_truefalse)
n_replace  <- sum(mask_truefalse, na.rm = TRUE)

cat("Celdas a reemplazar (TRUE en máscara):", n_replace, "de", n_total, 
    sprintf("(%.4f%%)\n", 100 * n_replace / n_total))

# Nota: si hay NA en incert_Fire51 en posiciones TRUE, se copiarán tal cual (NA).
incert_Fire60[mask_truefalse] <- incert_Fire51[mask_truefalse]

image.plot(lon, lat, incert_Fire60[,,2])
image.plot(lon, lat, incert_Fire51[,,2])

# --- 6) Guardar resultado (archivo nuevo para conservar el original) ---
out_replaced <- file.path(output_dir_RData, "BA_Incertidumbre_MRBA60.RData")
save(incert_Fire60, file = out_replaced)

cat("Guardado:\n  ", out_replaced, "\n")

# (Opcional) Si prefieres machacar el original, descomenta:
# save(incert_Fire60, file = incert60_file)

#Mascara final para producto solo grid donde FireCCI51 se mantiene ==1



## Visualización de comprobación
image.plot(mask_truefalse[, , 4])

##==============================================================================
## FINAL HARMONISATION STATUS MASK
##  NA = no MRBA60 burned-area value
##   0 = harmonisation applied
##   1 = FireCCI51 retained
##==============================================================================

load(file.path(output_dir_RData, "BA_MRBA60.RData"))

# stopifnot(identical(dim(Mask_union_global), dim(BA_FIRE60)))
dev.off()
## Keep original BA object untouched
BA_MRBA60_for_mask <- BA_FIRE60

## Treat BA = 0 as absence for this status mask
BA_MRBA60_for_mask[BA_MRBA60_for_mask == 0] <- NA

## Pixels with valid MRBA60 burned area
valid_MRBA60_BA <- !is.na(BA_MRBA60_for_mask) & BA_MRBA60_for_mask > 0

## Original mask:
## TRUE = FireCCI51 retained according to the original no-harmonisation mask
mask_FireCCI51_retained_raw <- Mask_union_global != 0

## Initialise final status mask as NA everywhere
Mask_MRBA60_harmonisation_status <- array(
  NA_integer_,
  dim = dim(Mask_union_global),
  dimnames = dimnames(Mask_union_global)
)

## 0 = harmonisation applied
Mask_MRBA60_harmonisation_status[valid_MRBA60_BA] <- 0L

## 1 = FireCCI51 retained,
## but ONLY where MRBA60 has actual burned area
Mask_MRBA60_harmonisation_status[
  valid_MRBA60_BA & mask_FireCCI51_retained_raw
] <- 1L

##==============================================================================
## CHECKS
##==============================================================================

## There should be no 0/1 where MRBA60 BA is NA
sum(!is.na(Mask_MRBA60_harmonisation_status) & is.na(BA_MRBA60_for_mask))

## Count classes
table(
  as.vector(Mask_MRBA60_harmonisation_status),
  useNA = "ifany"
)

## Check one month
table(
  as.vector(Mask_MRBA60_harmonisation_status[, , 1]),
  useNA = "ifany"
)


##==============================================================================
## QUICK PLOT
##==============================================================================

# dim(Mask_MRBA60_harmonisation_status)
# xx<-Mask_MRBA60_harmonisation_status[,,264]
# dim(xx)
# image.plot(
#   xx,
#   zlim = c(0, 1),
#   col = c("grey80", "black")
# )
# rm(xx);gc()
image.plot(
  Mask_MRBA60_harmonisation_status[, , 1],
  zlim = c(0, 1),
  col = c("grey80", "black")
)
dev.off()
##==============================================================================
## COUNT 0 / 1 AND PERCENTAGE OF VALUE 1
## NA values are excluded from the total
##==============================================================================
xx<-Mask_MRBA60_harmonisation_status[,,1:192]
x <- as.vector(xx)

n_0 <- sum(x == 0, na.rm = TRUE)
n_1 <- sum(x == 1, na.rm = TRUE)

n_total <- n_0 + n_1

pct_0 <- ifelse(n_total > 0, 100 * n_0 / n_total, NA_real_)
pct_1 <- ifelse(n_total > 0, 100 * n_1 / n_total, NA_real_)

mask_summary <- data.frame(
  n_0 = n_0,
  n_1 = n_1,
  n_total_valid = n_total,
  pct_0 = pct_0,
  pct_1 = pct_1
)

print(mask_summary)

##==============================================================================
## SAVE
##==============================================================================

rm(Mask_MRBA60_harmonisation_status);gc()
Mask_MRBA60_harmonisation_status<-xx
out_mask_final <- file.path(
  output_dir_RData,
  "MASK_MRBA60_HARMONISATION_STATUS.RData"
)

save(
  Mask_MRBA60_harmonisation_status,
  file = out_mask_final
)
