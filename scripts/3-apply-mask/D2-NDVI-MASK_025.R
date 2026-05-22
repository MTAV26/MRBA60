# ===============================
# CONFIGURACIÓN PREVIA
# ===============================

rm(list = ls())
graphics.off()
gc()

# Cargar librerías necesarias
library(ncdf4)  # leer archivos .nc
library(RANN)   # aplicar KNN con nn2()

# ============================================================================
# DIRECTORIOS ACTUALES
# ============================================================================

dir_ndvi <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/"
dir_adj  <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
dir_mask <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/"

# ============================================================================
# ARCHIVOS
# ============================================================================

file_ndvi <- paste0(dir_ndvi, "MOD13C2_NDVI_2003_2024_025.nc")

file_landsea <- paste0(
  dir_adj,
  "land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc"
)

file_firemask <- paste0(dir_mask, "FireMask_AF3030F.RData")

file_out <- paste0(
  dir_adj,
  "NDVI-2003_2024-MONTHLY-025-mask-landsea-KNN.RData"
)

# ============================================================================
# PARÁMETROS KNN
# ============================================================================

k_nn <- 8

# Distancia máxima permitida para buscar vecinos.
# El grid es de 0.25 grados.
# Por ejemplo:
#   4 celdas  = 1 grado aprox.
#   8 celdas  = 2 grados aprox.
#   12 celdas = 3 grados aprox.
max_dist_cells <- 8

# ============================================================================
# CARGAR LONGITUD Y LATITUD
# ============================================================================

load(paste0(dir_adj, "longitude.RData"))
load(paste0(dir_adj, "latitude.RData"))

# ============================================================================
# FECHAS: PERÍODO COMPLETO ACTUAL 2003-2024
# ============================================================================

dates_full <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")

anni <- 2003:2024
mesi <- rep(1:12, length(anni))
fechas <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")

inicio <- which(fechas == as.Date("2003-01-01"))
fin    <- which(fechas == as.Date("2024-12-01"))

# ============================================================================
# CARGAR NDVI DESDE NETCDF
# ============================================================================

nc <- nc_open(file_ndvi)

print(names(nc$var))

var_ndvi <- names(nc$var)[grepl("NDVI", names(nc$var), ignore.case = TRUE)][1]

if (is.na(var_ndvi)) {
  var_ndvi <- names(nc$var)[1]
}

cat("Variable NDVI usada:", var_ndvi, "\n")

NDVI_tot <- ncvar_get(nc, var_ndvi)

nc_close(nc)

cat("Dimensiones NDVI_tot:\n")
print(dim(NDVI_tot))

# ============================================================================
# COMPROBAR DIMENSIONES DEL NDVI
# ============================================================================

nrows <- dim(NDVI_tot)[1]
ncols <- dim(NDVI_tot)[2]
time_common <- dim(NDVI_tot)[3]

cat("Número de meses esperado:", length(dates_full), "\n")
cat("Número de meses en NDVI :", time_common, "\n")

if (time_common != length(dates_full)) {
  warning("El número de meses del NDVI no coincide con 2003-2024.")
}

# ============================================================================
# CARGAR MÁSCARA DE FUEGO
# ============================================================================

load(file_firemask)

status_matrix2_tot <- FireMask_AF3030F

cat("Dimensiones FireMask_AF3030F:\n")
print(dim(status_matrix2_tot))

cat("Valores FireMask_AF3030F:\n")
print(table(as.vector(status_matrix2_tot), useNA = "ifany"))

# ============================================================================
# CARGAR LAND/SEA MASK AJUSTADA
# ============================================================================

nc_land <- nc_open(file_landsea)

print(names(nc_land$var))

lon_land <- ncvar_get(nc_land, "lon")
lat_land <- ncvar_get(nc_land, "lat")

land_sea_mask <- ncvar_get(nc_land, "sftlf")

nc_close(nc_land)

cat("Dimensiones land_sea_mask:\n")
print(dim(land_sea_mask))

cat("Valores land_sea_mask:\n")
print(table(as.vector(land_sea_mask), useNA = "ifany"))


# ============================================================================
# COMPROBAR MÁSCARA DE PÍXELES A RELLENAR
# ============================================================================

# Solo se rellenan píxeles:
#   1) con fuego detectado
#   2) con NDVI ausente
#   3) clasificados como tierra

n_target_initial <- 0
n_target_ocean <- 0

for (tt in seq_len(time_common)) {
  
  n_target_initial <- n_target_initial +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(NDVI_tot[, , tt]) &
        land_sea_mask == 1,
      na.rm = TRUE
    )
  
  n_target_ocean <- n_target_ocean +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(NDVI_tot[, , tt]) &
        land_sea_mask == 0,
      na.rm = TRUE
    )
}

cat("Píxeles-tiempo a rellenar sobre tierra:", n_target_initial, "\n")
cat("Píxeles-tiempo con fuego y NDVI NA sobre océano:", n_target_ocean, "\n")

# Visualización rápida de un mes
mask_mes_1 <- status_matrix2_tot[, , 1] == 1 &
  is.na(NDVI_tot[, , 1]) &
  land_sea_mask == 1

image(mask_mes_1, main = "Píxeles a rellenar - mes 1")

# ============================================================================
# KNN PARA RELLENAR NDVI SOLO EN TIERRA Y CON DISTANCIA MÁXIMA
# ============================================================================

n_rellenos_total <- 0
n_sin_vecinos_total <- 0
n_fuera_distancia_total <- 0

fill_log <- data.frame(
  date = dates_full,
  n_missing_land = NA_integer_,
  n_valid_land = NA_integer_,
  n_filled = NA_integer_,
  n_no_donor = NA_integer_,
  n_far = NA_integer_,
  n_remaining = NA_integer_
)

for (tt in seq_len(time_common)) {
  
  cat("\nProcesando:", as.character(dates_full[tt]), "\n")
  
  ndvi_mes <- NDVI_tot[, , tt]
  fire_mes <- status_matrix2_tot[, , tt]
  
  # --------------------------------------------------------------------------
  # 1. Coordenadas válidas:
  #    NDVI no NA y tierra según land_sea_mask
  # --------------------------------------------------------------------------
  
  valid_xy <- which(
    !is.na(ndvi_mes) &
      land_sea_mask == 1,
    arr.ind = TRUE
  )
  
  # --------------------------------------------------------------------------
  # 2. Coordenadas a rellenar:
  #    fuego == 1, NDVI NA y tierra
  # --------------------------------------------------------------------------
  
  missing_xy <- which(
    fire_mes == 1 &
      is.na(ndvi_mes) &
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
    fill_log$n_far[tt] <- 0
    fill_log$n_remaining[tt] <- 0
    
    next
  }
  
  if (n_valid == 0) {
    
    warning("No hay píxeles donantes válidos sobre tierra para ", dates_full[tt])
    
    fill_log$n_missing_land[tt] <- n_missing
    fill_log$n_valid_land[tt] <- 0
    fill_log$n_filled[tt] <- 0
    fill_log$n_no_donor[tt] <- n_missing
    fill_log$n_far[tt] <- 0
    fill_log$n_remaining[tt] <- n_missing
    
    n_sin_vecinos_total <- n_sin_vecinos_total + n_missing
    
    next
  }
  
  # Si hay menos donantes que k_nn, se ajusta k
  k_eff <- min(k_nn, n_valid)
  
  # --------------------------------------------------------------------------
  # 3. Buscar vecinos KNN dentro del mismo mes
  # --------------------------------------------------------------------------
  
  nn <- nn2(
    data  = valid_xy,
    query = missing_xy,
    k     = k_eff
  )
  
  n_rellenos_mes <- 0
  n_far_mes <- 0
  
  for (i in seq_len(nrow(missing_xy))) {
    
    r0 <- missing_xy[i, 1]
    c0 <- missing_xy[i, 2]
    
    neigh_idx <- nn$nn.idx[i, ]
    neigh_dist <- nn$nn.dists[i, ]
    
    # ------------------------------------------------------------------------
    # 4. Mantener solo vecinos dentro de la distancia máxima
    # ------------------------------------------------------------------------
    
    neigh_ok <- neigh_dist <= max_dist_cells
    
    if (!any(neigh_ok)) {
      n_far_mes <- n_far_mes + 1
      next
    }
    
    neigh_xy <- valid_xy[neigh_idx[neigh_ok], , drop = FALSE]
    
    vals <- apply(
      neigh_xy,
      1,
      function(x) ndvi_mes[x[1], x[2]]
    )
    
    vals_validos <- vals[!is.na(vals)]
    
    if (length(vals_validos) > 0) {
      ndvi_mes[r0, c0] <- mean(vals_validos)
      n_rellenos_mes <- n_rellenos_mes + 1
    }
  }
  
  # Actualizar el cubo NDVI
  NDVI_tot[, , tt] <- ndvi_mes
  
  # Comprobación mensual
  n_remaining_mes <- sum(
    fire_mes == 1 &
      is.na(NDVI_tot[, , tt]) &
      land_sea_mask == 1,
    na.rm = TRUE
  )
  
  n_rellenos_total <- n_rellenos_total + n_rellenos_mes
  n_fuera_distancia_total <- n_fuera_distancia_total + n_far_mes
  
  fill_log$n_missing_land[tt] <- n_missing
  fill_log$n_valid_land[tt] <- n_valid
  fill_log$n_filled[tt] <- n_rellenos_mes
  fill_log$n_no_donor[tt] <- 0
  fill_log$n_far[tt] <- n_far_mes
  fill_log$n_remaining[tt] <- n_remaining_mes
  
  cat("Puntos rellenados:", n_rellenos_mes, "\n")
  cat("Puntos sin vecinos dentro de distancia máxima:", n_far_mes, "\n")
  cat("Puntos restantes:", n_remaining_mes, "\n")
  
  rm(
    ndvi_mes,
    fire_mes,
    valid_xy,
    missing_xy,
    nn
  )
  
  gc()
}

# ============================================================================
# COMPROBACIÓN FINAL
# ============================================================================

n_remaining_final <- 0

for (tt in seq_len(time_common)) {
  
  n_remaining_final <- n_remaining_final +
    sum(
      status_matrix2_tot[, , tt] == 1 &
        is.na(NDVI_tot[, , tt]) &
        land_sea_mask == 1,
      na.rm = TRUE
    )
}

cat("\n================ RESUMEN FINAL ================\n")
cat("Píxeles-tiempo iniciales a rellenar sobre tierra:", n_target_initial, "\n")
cat("Píxeles-tiempo rellenados:", n_rellenos_total, "\n")
cat("Píxeles-tiempo sin vecinos próximos:", n_fuera_distancia_total, "\n")
cat("Píxeles-tiempo restantes sobre tierra:", n_remaining_final, "\n")

cat("\nResumen mensual:\n")
print(fill_log)

# ============================================================================
# GUARDAR RESULTADO FINAL
# ============================================================================

save(
  lon,
  lat,
  NDVI_tot,
  file = file_out
)

cat("\nArchivo guardado en:\n")
cat(file_out, "\n")
# ============================================================================
# COMPARACIÓN VISUAL Y RESUMEN
# ============================================================================

NDVI_new <- NDVI_tot
rm(NDVI_tot)

image(NDVI_new[, , 1], main = "NDVI rellenado - mes 1")
summary(as.vector(NDVI_new))

cat("\nArchivo guardado en:\n")
cat(file_out, "\n")
