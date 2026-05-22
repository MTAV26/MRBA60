# ============================================================
# VPD 2003-2024: RELLENO KNN SOBRE TIERRA + PLOTS DE CONTROL
# ============================================================

rm(list = ls())
graphics.off()
gc()

# ============================================================
# PAQUETES NECESARIOS
# ============================================================

library(ncdf4)
library(RANN)

library(ggplot2)
library(viridis)
library(patchwork)
library(rnaturalearth)
library(sf)
library(grid)

# ============================================================
# DIRECTORIOS
# ============================================================

dir_vpd  <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/"
dir_adj  <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
dir_mask <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/"

dir_plot <- paste0(dir_adj, "plots/")

if (!dir.exists(dir_plot)) {
  dir.create(dir_plot, recursive = TRUE)
}

# ============================================================
# ARCHIVOS
# ============================================================

file_vpd <- paste0(dir_vpd, "vpd_2003-2024_0.25deg_bil.nc")

file_landsea <- paste0(
  dir_adj,
  "land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc"
)

file_firemask <- paste0(dir_mask, "FireMask_AF3030F.RData")

file_out <- paste0(
  dir_adj,
  "VPD-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"
)

file_out_with_original <- paste0(
  dir_adj,
  "VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-with-original.RData"
)

file_log_csv <- paste0(
  dir_adj,
  "VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-fill_log.csv"
)

file_summary_rdata <- paste0(
  dir_adj,
  "VPD-2003_2024-MONTHLY-025-mask-landsea-KNN-summary.RData"
)

# ============================================================
# PARÁMETROS KNN
# ============================================================

k_nn <- 8

# Radio máximo inicial.
# Grid de 0.25º:
#   8 celdas = 2º aprox.
max_dist_cells <- 8
grid_res_deg <- 0.25

# ============================================================
# CARGAR LONGITUD Y LATITUD
# ============================================================

load(paste0(dir_adj, "longitude.RData"))
load(paste0(dir_adj, "latitude.RData"))

# Normalizar nombres
if (exists("longitude")) lon <- longitude
if (exists("latitude"))  lat <- latitude

if (!exists("lon")) stop("No se encontró objeto lon o longitude.")
if (!exists("lat")) stop("No se encontró objeto lat o latitude.")

cat("Longitudes:", length(lon), "\n")
cat("Latitudes :", length(lat), "\n")

# ============================================================
# FECHAS
# ============================================================

dates_full <- seq(
  as.Date("2003-01-01"),
  as.Date("2024-12-01"),
  by = "month"
)

anni <- 2003:2024
mesi <- rep(1:12, length(anni))
fechas <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")

inicio <- which(fechas == as.Date("2003-01-01"))
fin    <- which(fechas == as.Date("2024-12-01"))

# ============================================================
# CARGAR VPD DESDE NETCDF
# ============================================================

nc <- nc_open(file_vpd)

cat("\nVariables disponibles en VPD NetCDF:\n")
print(names(nc$var))

var_vpd <- names(nc$var)[grepl("VPD|vpd", names(nc$var), ignore.case = TRUE)][1]

if (is.na(var_vpd)) {
  var_vpd <- names(nc$var)[1]
}

cat("Variable VPD usada:", var_vpd, "\n")

VPD_tot <- ncvar_get(nc, var_vpd)

nc_close(nc)

cat("Dimensiones VPD_tot:\n")
print(dim(VPD_tot))

# Guardar copia original antes del relleno KNN
VPD_original <- VPD_tot

# ============================================================
# COMPROBAR DIMENSIONES DEL VPD
# ============================================================

nrows <- dim(VPD_tot)[1]
ncols <- dim(VPD_tot)[2]
time_common <- dim(VPD_tot)[3]

cat("Número de meses esperado:", length(dates_full), "\n")
cat("Número de meses en VPD  :", time_common, "\n")

if (time_common != length(dates_full)) {
  warning("El número de meses del VPD no coincide con 2003-2024.")
}

# ============================================================
# CARGAR MÁSCARA DE FUEGO
# ============================================================

load(file_firemask)

status_matrix2_tot <- FireMask_AF3030F

cat("\nDimensiones FireMask_AF3030F:\n")
print(dim(status_matrix2_tot))

cat("Valores FireMask_AF3030F:\n")
print(table(as.vector(status_matrix2_tot), useNA = "ifany"))

# ============================================================
# CARGAR LAND/SEA MASK AJUSTADA
# ============================================================

nc_land <- nc_open(file_landsea)

cat("\nVariables disponibles en land/sea NetCDF:\n")
print(names(nc_land$var))

lon_land <- ncvar_get(nc_land, "lon")
lat_land <- ncvar_get(nc_land, "lat")

land_sea_mask <- ncvar_get(nc_land, "sftlf")

nc_close(nc_land)

cat("Dimensiones land_sea_mask:\n")
print(dim(land_sea_mask))

cat("Valores land_sea_mask:\n")
print(table(as.vector(land_sea_mask), useNA = "ifany"))

# ============================================================
# COMPROBAR DIMENSIONES
# ============================================================

if (!all(dim(VPD_tot)[1:2] == dim(status_matrix2_tot)[1:2])) {
  stop(
    "Las dimensiones espaciales de VPD_tot y FireMask_AF3030F no coinciden.\n",
    "VPD_tot: ", paste(dim(VPD_tot)[1:2], collapse = " x "), "\n",
    "FireMask: ", paste(dim(status_matrix2_tot)[1:2], collapse = " x ")
  )
}

if (!all(dim(VPD_tot)[1:2] == dim(land_sea_mask))) {
  stop(
    "Las dimensiones espaciales de VPD_tot y land_sea_mask no coinciden.\n",
    "VPD_tot: ", paste(dim(VPD_tot)[1:2], collapse = " x "), "\n",
    "land_sea_mask: ", paste(dim(land_sea_mask), collapse = " x ")
  )
}

if (dim(VPD_tot)[3] != dim(status_matrix2_tot)[3]) {
  stop(
    "La dimensión temporal de VPD_tot y FireMask_AF3030F no coincide.\n",
    "VPD_tot time: ", dim(VPD_tot)[3], "\n",
    "FireMask time: ", dim(status_matrix2_tot)[3]
  )
}

# ============================================================
# COMPROBAR PÍXELES A RELLENAR
# ============================================================

n_target_initial <- 0
n_target_ocean <- 0

for (tt in seq_len(time_common)) {
  
  n_target_initial <- n_target_initial +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(VPD_tot[, , tt]) &
        land_sea_mask == 1,
      na.rm = TRUE
    )
  
  n_target_ocean <- n_target_ocean +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(VPD_tot[, , tt]) &
        land_sea_mask == 0,
      na.rm = TRUE
    )
}

cat("\nPíxeles-tiempo a rellenar sobre tierra:", n_target_initial, "\n")
cat("Píxeles-tiempo con fuego y VPD NA sobre océano:", n_target_ocean, "\n")

mask_mes_1 <- status_matrix2_tot[, , 1] == 1 &
  is.na(VPD_tot[, , 1]) &
  land_sea_mask == 1

image(mask_mes_1, main = "Píxeles VPD a rellenar - mes 1")

# ============================================================
# KNN PARA RELLENAR VPD SOLO EN TIERRA
# Primer intento: vecinos dentro de max_dist_cells
# Respaldo: si no hay vecinos dentro del radio, usar los k vecinos
# terrestres válidos más cercanos disponibles y registrar distancia.
# ============================================================

n_rellenos_total <- 0
n_sin_vecinos_total <- 0
n_rellenos_dentro_total <- 0
n_rellenos_fuera_total <- 0

distancias_fuera_total_min <- c()
distancias_fuera_total_mean <- c()
distancias_fuera_total_max <- c()

fill_log <- data.frame(
  date = dates_full,
  n_missing_land = NA_integer_,
  n_valid_land = NA_integer_,
  n_filled = NA_integer_,
  n_no_donor = NA_integer_,
  n_filled_within_maxdist = NA_integer_,
  n_filled_beyond_maxdist = NA_integer_,
  min_dist_beyond_cells = NA_real_,
  mean_dist_beyond_cells = NA_real_,
  max_dist_beyond_cells = NA_real_,
  min_dist_beyond_degrees = NA_real_,
  mean_dist_beyond_degrees = NA_real_,
  max_dist_beyond_degrees = NA_real_,
  n_remaining = NA_integer_
)

for (tt in seq_len(time_common)) {
  
  cat("\nProcesando:", as.character(dates_full[tt]), "\n")
  
  vpd_mes <- VPD_tot[, , tt]
  fire_mes <- status_matrix2_tot[, , tt]
  
  # Donantes: VPD válido y tierra
  valid_xy <- which(
    !is.na(vpd_mes) &
      land_sea_mask == 1,
    arr.ind = TRUE
  )
  
  # Objetivo: fuego, VPD NA y tierra
  missing_xy <- which(
    fire_mes == 1 &
      is.na(vpd_mes) &
      land_sea_mask == 1,
    arr.ind = TRUE
  )
  
  n_missing <- nrow(missing_xy)
  n_valid <- nrow(valid_xy)
  
  cat("Puntos a rellenar sobre tierra:", n_missing, "\n")
  cat("Puntos válidos disponibles sobre tierra:", n_valid, "\n")
  
  if (n_missing == 0) {
    
    fill_log$n_missing_land[tt] <- 0
    fill_log$n_valid_land[tt] <- n_valid
    fill_log$n_filled[tt] <- 0
    fill_log$n_no_donor[tt] <- 0
    fill_log$n_filled_within_maxdist[tt] <- 0
    fill_log$n_filled_beyond_maxdist[tt] <- 0
    fill_log$n_remaining[tt] <- 0
    
    next
  }
  
  if (n_valid == 0) {
    
    warning("No hay píxeles donantes válidos sobre tierra para ", dates_full[tt])
    
    fill_log$n_missing_land[tt] <- n_missing
    fill_log$n_valid_land[tt] <- 0
    fill_log$n_filled[tt] <- 0
    fill_log$n_no_donor[tt] <- n_missing
    fill_log$n_filled_within_maxdist[tt] <- 0
    fill_log$n_filled_beyond_maxdist[tt] <- 0
    fill_log$n_remaining[tt] <- n_missing
    
    n_sin_vecinos_total <- n_sin_vecinos_total + n_missing
    
    next
  }
  
  k_eff <- min(k_nn, n_valid)
  
  nn <- nn2(
    data  = valid_xy,
    query = missing_xy,
    k     = k_eff
  )
  
  n_rellenos_mes <- 0
  n_far_mes <- 0
  
  distancias_fuera_mes_min <- c()
  distancias_fuera_mes_mean <- c()
  distancias_fuera_mes_max <- c()
  
  for (i in seq_len(nrow(missing_xy))) {
    
    r0 <- missing_xy[i, 1]
    c0 <- missing_xy[i, 2]
    
    neigh_idx <- nn$nn.idx[i, ]
    neigh_dist <- nn$nn.dists[i, ]
    
    # Primero vecinos dentro del radio máximo
    neigh_ok <- neigh_dist <= max_dist_cells
    
    relleno_fuera_distancia <- FALSE
    distancia_min_necesaria <- NA_real_
    distancia_mean_usada <- NA_real_
    distancia_max_usada <- NA_real_
    
    # Si no hay vecinos dentro del radio,
    # usar los k vecinos terrestres válidos más cercanos disponibles.
    if (!any(neigh_ok)) {
      
      relleno_fuera_distancia <- TRUE
      
      neigh_ok <- rep(TRUE, length(neigh_dist))
      
      distancia_min_necesaria <- min(neigh_dist, na.rm = TRUE)
      distancia_mean_usada <- mean(neigh_dist, na.rm = TRUE)
      distancia_max_usada <- max(neigh_dist, na.rm = TRUE)
    }
    
    neigh_xy <- valid_xy[neigh_idx[neigh_ok], , drop = FALSE]
    
    vals <- apply(
      neigh_xy,
      1,
      function(x) vpd_mes[x[1], x[2]]
    )
    
    vals_validos <- vals[!is.na(vals)]
    
    if (length(vals_validos) > 0) {
      
      vpd_mes[r0, c0] <- mean(vals_validos)
      n_rellenos_mes <- n_rellenos_mes + 1
      
      if (relleno_fuera_distancia) {
        
        n_far_mes <- n_far_mes + 1
        
        distancias_fuera_mes_min <- c(
          distancias_fuera_mes_min,
          distancia_min_necesaria
        )
        
        distancias_fuera_mes_mean <- c(
          distancias_fuera_mes_mean,
          distancia_mean_usada
        )
        
        distancias_fuera_mes_max <- c(
          distancias_fuera_mes_max,
          distancia_max_usada
        )
      }
    }
  }
  
  # Actualizar cubo VPD
  VPD_tot[, , tt] <- vpd_mes
  
  n_remaining_mes <- sum(
    fire_mes == 1 &
      is.na(VPD_tot[, , tt]) &
      land_sea_mask == 1,
    na.rm = TRUE
  )
  
  n_rellenos_dentro_mes <- n_rellenos_mes - n_far_mes
  
  n_rellenos_total <- n_rellenos_total + n_rellenos_mes
  n_rellenos_dentro_total <- n_rellenos_dentro_total + n_rellenos_dentro_mes
  n_rellenos_fuera_total <- n_rellenos_fuera_total + n_far_mes
  
  distancias_fuera_total_min <- c(
    distancias_fuera_total_min,
    distancias_fuera_mes_min
  )
  
  distancias_fuera_total_mean <- c(
    distancias_fuera_total_mean,
    distancias_fuera_mes_mean
  )
  
  distancias_fuera_total_max <- c(
    distancias_fuera_total_max,
    distancias_fuera_mes_max
  )
  
  fill_log$n_missing_land[tt] <- n_missing
  fill_log$n_valid_land[tt] <- n_valid
  fill_log$n_filled[tt] <- n_rellenos_mes
  fill_log$n_no_donor[tt] <- 0
  fill_log$n_filled_within_maxdist[tt] <- n_rellenos_dentro_mes
  fill_log$n_filled_beyond_maxdist[tt] <- n_far_mes
  
  if (length(distancias_fuera_mes_min) > 0) {
    
    fill_log$min_dist_beyond_cells[tt] <- min(distancias_fuera_mes_min, na.rm = TRUE)
    fill_log$mean_dist_beyond_cells[tt] <- mean(distancias_fuera_mes_mean, na.rm = TRUE)
    fill_log$max_dist_beyond_cells[tt] <- max(distancias_fuera_mes_max, na.rm = TRUE)
    
    fill_log$min_dist_beyond_degrees[tt] <- min(distancias_fuera_mes_min, na.rm = TRUE) * grid_res_deg
    fill_log$mean_dist_beyond_degrees[tt] <- mean(distancias_fuera_mes_mean, na.rm = TRUE) * grid_res_deg
    fill_log$max_dist_beyond_degrees[tt] <- max(distancias_fuera_mes_max, na.rm = TRUE) * grid_res_deg
  }
  
  fill_log$n_remaining[tt] <- n_remaining_mes
  
  cat("Puntos rellenados:", n_rellenos_mes, "\n")
  cat("Puntos rellenados dentro de distancia máxima:", n_rellenos_dentro_mes, "\n")
  cat("Puntos rellenados fuera de distancia máxima:", n_far_mes, "\n")
  
  if (length(distancias_fuera_mes_min) > 0) {
    
    cat("Distancia mínima necesaria fuera del radio, en celdas:",
        min(distancias_fuera_mes_min, na.rm = TRUE), "\n")
    
    cat("Distancia media usada fuera del radio, en celdas:",
        mean(distancias_fuera_mes_mean, na.rm = TRUE), "\n")
    
    cat("Distancia máxima usada fuera del radio, en celdas:",
        max(distancias_fuera_mes_max, na.rm = TRUE), "\n")
    
    cat("Distancia máxima usada fuera del radio, en grados aprox.:",
        max(distancias_fuera_mes_max, na.rm = TRUE) * grid_res_deg, "\n")
  }
  
  cat("Puntos restantes:", n_remaining_mes, "\n")
  
  rm(vpd_mes, fire_mes, valid_xy, missing_xy, nn)
  gc()
}

# ============================================================
# COMPROBACIÓN FINAL
# ============================================================

n_remaining_final <- 0

for (tt in seq_len(time_common)) {
  
  n_remaining_final <- n_remaining_final +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(VPD_tot[, , tt]) &
        land_sea_mask == 1,
      na.rm = TRUE
    )
}

cat("\n================ RESUMEN FINAL ================\n")
cat("Píxeles-tiempo iniciales a rellenar sobre tierra:", n_target_initial, "\n")
cat("Píxeles-tiempo rellenados:", n_rellenos_total, "\n")
cat("Píxeles-tiempo rellenados dentro de distancia máxima:", n_rellenos_dentro_total, "\n")
cat("Píxeles-tiempo rellenados fuera de distancia máxima:", n_rellenos_fuera_total, "\n")
cat("Píxeles-tiempo sin donantes terrestres válidos:", n_sin_vecinos_total, "\n")
cat("Píxeles-tiempo restantes sobre tierra:", n_remaining_final, "\n")

if (length(distancias_fuera_total_min) > 0) {
  
  cat("\nDistancias necesarias para los casos fuera del radio máximo:\n")
  
  cat("Distancia mínima al vecino más cercano, en celdas:",
      min(distancias_fuera_total_min, na.rm = TRUE), "\n")
  
  cat("Distancia media de los vecinos usados, en celdas:",
      mean(distancias_fuera_total_mean, na.rm = TRUE), "\n")
  
  cat("Distancia máxima de los vecinos usados, en celdas:",
      max(distancias_fuera_total_max, na.rm = TRUE), "\n")
  
  cat("Distancia mínima al vecino más cercano, en grados aprox.:",
      min(distancias_fuera_total_min, na.rm = TRUE) * grid_res_deg, "\n")
  
  cat("Distancia media de los vecinos usados, en grados aprox.:",
      mean(distancias_fuera_total_mean, na.rm = TRUE) * grid_res_deg, "\n")
  
  cat("Distancia máxima de los vecinos usados, en grados aprox.:",
      max(distancias_fuera_total_max, na.rm = TRUE) * grid_res_deg, "\n")
}

cat("\nResumen mensual:\n")
print(fill_log)

# ============================================================
# GUARDAR LOG Y RESUMEN
# ============================================================

write.csv(
  fill_log,
  file = file_log_csv,
  row.names = FALSE
)

summary_knn <- list(
  variable = "VPD",
  period = "2003-2024",
  k_nn = k_nn,
  max_dist_cells = max_dist_cells,
  grid_res_deg = grid_res_deg,
  max_dist_degrees = max_dist_cells * grid_res_deg,
  n_target_initial = n_target_initial,
  n_target_ocean = n_target_ocean,
  n_rellenos_total = n_rellenos_total,
  n_rellenos_dentro_total = n_rellenos_dentro_total,
  n_rellenos_fuera_total = n_rellenos_fuera_total,
  n_sin_vecinos_total = n_sin_vecinos_total,
  n_remaining_final = n_remaining_final,
  dist_min_beyond_cells = ifelse(
    length(distancias_fuera_total_min) > 0,
    min(distancias_fuera_total_min, na.rm = TRUE),
    NA
  ),
  dist_mean_beyond_cells = ifelse(
    length(distancias_fuera_total_mean) > 0,
    mean(distancias_fuera_total_mean, na.rm = TRUE),
    NA
  ),
  dist_max_beyond_cells = ifelse(
    length(distancias_fuera_total_max) > 0,
    max(distancias_fuera_total_max, na.rm = TRUE),
    NA
  ),
  dist_min_beyond_degrees = ifelse(
    length(distancias_fuera_total_min) > 0,
    min(distancias_fuera_total_min, na.rm = TRUE) * grid_res_deg,
    NA
  ),
  dist_mean_beyond_degrees = ifelse(
    length(distancias_fuera_total_mean) > 0,
    mean(distancias_fuera_total_mean, na.rm = TRUE) * grid_res_deg,
    NA
  ),
  dist_max_beyond_degrees = ifelse(
    length(distancias_fuera_total_max) > 0,
    max(distancias_fuera_total_max, na.rm = TRUE) * grid_res_deg,
    NA
  )
)

save(
  summary_knn,
  fill_log,
  file = file_summary_rdata
)

cat("\nLog mensual guardado en:\n")
cat(file_log_csv, "\n")

cat("\nResumen KNN guardado en:\n")
cat(file_summary_rdata, "\n")

# ============================================================
# GUARDAR RESULTADOS
# ============================================================

# RData principal: solo lon, lat y VPD_tot
save(
  lon,
  lat,
  VPD_tot,
  file = file_out
)

# RData auxiliar con original para reproducir plots
save(
  lon,
  lat,
  VPD_original,
  VPD_tot,
  status_matrix2_tot,
  land_sea_mask,
  dates_full,
  fill_log,
  summary_knn,
  file = file_out_with_original
)

cat("\nArchivo principal guardado en:\n")
cat(file_out, "\n")

cat("\nArchivo auxiliar con VPD original guardado en:\n")
cat(file_out_with_original, "\n")

# ============================================================
# PLOTS DE CONTROL
# ============================================================

cat("\nResumen VPD original:\n")
print(summary(as.vector(VPD_original)))

cat("\nResumen VPD rellenado:\n")
print(summary(as.vector(VPD_tot)))

cat("\nResumen FireMask:\n")
print(table(as.vector(status_matrix2_tot), useNA = "ifany"))

# Mes a representar
tt_plot <- 1
plot_date <- dates_full[tt_plot]

cat("\nMes representado:\n")
print(plot_date)

# Región: SE Asia / Indonesia / Filipinas
xrange <- c(92, 130)
yrange <- c(-7, 17)

# ============================================================
# DATA FRAMES
# ============================================================

df_vpd_original <- data.frame(
  expand.grid(lon = lon, lat = lat),
  VPD = as.vector(VPD_original[, , tt_plot])
)

df_vpd_filled <- data.frame(
  expand.grid(lon = lon, lat = lat),
  VPD = as.vector(VPD_tot[, , tt_plot])
)

df_mask <- data.frame(
  expand.grid(lon = lon, lat = lat),
  mask = as.vector(status_matrix2_tot[, , tt_plot])
)

filled_cells <- is.na(VPD_original[, , tt_plot]) &
  !is.na(VPD_tot[, , tt_plot]) &
  status_matrix2_tot[, , tt_plot] == 1 &
  land_sea_mask == 1

df_filled_cells <- data.frame(
  expand.grid(lon = lon, lat = lat),
  filled = as.integer(as.vector(filled_cells))
)

# ============================================================
# FRONTERAS
# ============================================================

world <- ne_countries(scale = "medium", returnclass = "sf")

region_bbox <- st_as_sfc(
  st_bbox(
    c(
      xmin = xrange[1],
      xmax = xrange[2],
      ymin = yrange[1],
      ymax = yrange[2]
    ),
    crs = st_crs(world)
  )
)

world_crop <- st_crop(world, region_bbox)

# ============================================================
# TEMA Y ESCALA COMÚN
# ============================================================

border_theme <- theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  plot.margin  = unit(c(5, 5, 5, 5), "pt"),
  plot.title   = element_text(hjust = 0.5, size = 15),
  axis.title   = element_text(size = 12),
  axis.text    = element_text(size = 10),
  legend.title = element_text(size = 11),
  legend.text  = element_text(size = 9)
)

vpd_min <- min(
  c(df_vpd_original$VPD, df_vpd_filled$VPD),
  na.rm = TRUE
)

vpd_max <- max(
  c(df_vpd_original$VPD, df_vpd_filled$VPD),
  na.rm = TRUE
)

# ============================================================
# PLOT 1: VPD ORIGINAL
# ============================================================

p_original <- ggplot() +
  
  geom_raster(
    data = subset(df_mask, mask == 1),
    aes(x = lon, y = lat),
    fill = "red"
  ) +
  
  geom_raster(
    data = df_vpd_original,
    aes(x = lon, y = lat, fill = VPD)
  ) +
  
  geom_sf(
    data = world_crop,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    na.value = "transparent",
    limits = c(vpd_min, vpd_max)
  ) +
  
  coord_sf(
    xlim = xrange,
    ylim = yrange,
    expand = FALSE
  ) +
  
  labs(
    title = paste0(
      "a) Original VPD – ",
      format(plot_date, "%B %Y"),
      "\n(fire mask below)"
    ),
    x = "Longitude",
    y = "Latitude",
    fill = "VPD"
  ) +
  
  theme_minimal() +
  border_theme

# ============================================================
# PLOT 2: VPD RELLENADO
# ============================================================

p_filled <- ggplot() +
  
  geom_raster(
    data = subset(df_mask, mask == 1),
    aes(x = lon, y = lat),
    fill = "red"
  ) +
  
  geom_raster(
    data = df_vpd_filled,
    aes(x = lon, y = lat, fill = VPD)
  ) +
  
  geom_sf(
    data = world_crop,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  
  scale_fill_viridis_c(
    option = "viridis",
    na.value = "transparent",
    limits = c(vpd_min, vpd_max)
  ) +
  
  coord_sf(
    xlim = xrange,
    ylim = yrange,
    expand = FALSE
  ) +
  
  labs(
    title = paste0(
      "b) VPD after k-NN filling – ",
      format(plot_date, "%B %Y"),
      "\n(fire mask below)"
    ),
    x = "Longitude",
    y = "Latitude",
    fill = "VPD"
  ) +
  
  theme_minimal() +
  border_theme

combined_plot <- p_original + p_filled

print(combined_plot)

file_plot_pdf <- paste0(
  dir_plot,
  "VPD_original_vs_KNN_filled_SEAsia_",
  format(plot_date, "%Y_%m"),
  ".pdf"
)

file_plot_png <- paste0(
  dir_plot,
  "VPD_original_vs_KNN_filled_SEAsia_",
  format(plot_date, "%Y_%m"),
  ".png"
)

ggsave(
  filename = file_plot_pdf,
  plot = combined_plot,
  width = 12,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file_plot_png,
  plot = combined_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)

cat("\nFigura comparativa PDF guardada en:\n")
cat(file_plot_pdf, "\n")

cat("\nFigura comparativa PNG guardada en:\n")
cat(file_plot_png, "\n")

# ============================================================
# PLOT 3: CELDAS RELLENADAS POR KNN
# ============================================================

p_filled_cells <- ggplot() +
  
  geom_raster(
    data = subset(df_filled_cells, filled == 1),
    aes(x = lon, y = lat),
    fill = "red"
  ) +
  
  geom_sf(
    data = world_crop,
    fill = NA,
    color = "black",
    linewidth = 0.3
  ) +
  
  coord_sf(
    xlim = xrange,
    ylim = yrange,
    expand = FALSE
  ) +
  
  labs(
    title = paste0(
      "VPD cells filled by k-NN – ",
      format(plot_date, "%B %Y")
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  
  theme_minimal() +
  border_theme

print(p_filled_cells)

file_plot_filled_cells_pdf <- paste0(
  dir_plot,
  "VPD_cells_filled_by_KNN_SEAsia_",
  format(plot_date, "%Y_%m"),
  ".pdf"
)

file_plot_filled_cells_png <- paste0(
  dir_plot,
  "VPD_cells_filled_by_KNN_SEAsia_",
  format(plot_date, "%Y_%m"),
  ".png"
)

ggsave(
  filename = file_plot_filled_cells_pdf,
  plot = p_filled_cells,
  width = 6,
  height = 6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file_plot_filled_cells_png,
  plot = p_filled_cells,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300
)

cat("\nMapa de celdas rellenadas PDF guardado en:\n")
cat(file_plot_filled_cells_pdf, "\n")

cat("\nMapa de celdas rellenadas PNG guardado en:\n")
cat(file_plot_filled_cells_png, "\n")

# ============================================================
# FIN
# ============================================================

cat("\nProceso terminado correctamente.\n")