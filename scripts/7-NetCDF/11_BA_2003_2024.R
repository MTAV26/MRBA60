

rm(list = ls())
graphics.off()
gc()

# output_dir_RData<-"C:/Users/migue/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-2/RData/"
Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
output_dir_RData<-"/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/"

#MASK TRUE/fALSE
load(paste0(output_dir_RData, "MASK_NOHARMONISED.RData"))
# DATOS BA
load(paste0(output_dir_RData, "BA_MRBA60.RData"))
# load("/mnt/disco6tb//Dropbox/UAH/FireCCI60/FireCCI51_2001_2022_0.25degree-download.RData")
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")

dim(Mask_union_global)
dim(BA_FIRE60)
Fire51=f51
rm(f51);gc()
dim(Fire51)



# 1) Asegurar la máscara como lógica y sin NA
Mask_union_global <- (Mask_union_global == TRUE) | (Mask_union_global == 1)
Mask_union_global[is.na(Mask_union_global)] <- FALSE

stopifnot(identical(dim(Mask_union_global), dim(BA_FIRE60)))


Fire51_2003_2022 <- Fire51 /1e6
# image.plot(Fire51_2003_2022[,,1])
# Fire51_2003_2022[Fire51_2003_2022==0]=NA
stopifnot(identical(dim(Fire51_2003_2022), dim(BA_FIRE60)))
# 3) Construir FIRECCI60_FINAL: donde mask==TRUE usamos Fire51; si mask==FALSE, mantenemos BA_FIRE60
FIRECCI60_FINAL <- BA_FIRE60
FIRECCI60_FINAL[Mask_union_global] <- Fire51_2003_2022[Mask_union_global]
summary(as.vector(FIRECCI60_FINAL))





rm(Fire51, Fire51_2003_2022); gc()
# ---------------------------
# COMPROBAR MINIMO Y MAXIMO POSIBLE
# ---------------------------

# file_mask  <- file.path(output_dir_RData, paste0("F51_maskAboveTope_GLOBAL_", Modelo, "_2003_2022.RData"))
load(paste0(output_dir_RData, "/F51_maskAboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData"))
load(paste0(output_dir_RData, "F51_aboveTope_GLOBAL_B1-MRBA60-2003-2024_2003_2024.RData"))

# Asegura lógicas
Mask_Fire51_aboveTope_global <- (Mask_Fire51_aboveTope_global == 1) | (Mask_Fire51_aboveTope_global == TRUE)
Mask_Fire51_aboveTope_global[is.na(Mask_Fire51_aboveTope_global)] <- FALSE

# 1) parte de una copia completa -> así 193:240 ya queda relleno
BA_qm_global_co <- FIRECCI60_FINAL   # si quieres, usa BA_harmonised si ese es el baseline correcto

# 2) aplica solo en 2003–2018 (1:192)
for (t in 1:192) {
  m <- Mask_Fire51_aboveTope_global[,,t]
  if (!any(m)) next  # no hay nada que cambiar, pero el slice ya está copiado
  BA_qm_global_co[,,t][m] <- BA_Fire51_aboveTope_global[,,t][m]
}

# dim(global_BA_FireHarmonized_full)
cell_area_constant <- (110.57 * 0.25) * (111.32 * 0.25)  # ≈ 769.29 km² en el ecuador
area_by_row <- cell_area_constant * cos(lat * pi/180)
nrow_grid <- length(lat)
ncol_grid <- 1440
area_matrix <- matrix(area_by_row, ncol = ncol_grid, nrow = nrow_grid, byrow = FALSE)
area_matrix <- t(area_matrix)
image(lon, lat, area_matrix)
dim(area_matrix)

# --- Comprobación de excedencias BA > área de celda ---
check_ba_exceeds_area <- function(ba_arr, area_mat, lon, lat, tol = 0) {
  stopifnot(length(dim(ba_arr)) == 3)
  d <- dim(ba_arr)  # c(nx, ny, nt)
  stopifnot(identical(dim(area_mat), d[1:2]))
  stopifnot(length(lon) == d[1], length(lat) == d[2])
  area_3d <- array(rep(area_mat, d[3]), dim = d)  # Expandir área a 3D para comparar con BA por tiempo
  # Máscara de excedencias (permitiendo tolerancia 'tol', p.ej. tol = 1e-6)
  exceed_mask <- (ba_arr > (area_3d + tol)) & !is.na(ba_arr) & !is.na(area_3d)
  any_exceed <- any(exceed_mask, na.rm = TRUE)
  total_cells_time <- sum(!is.na(ba_arr))
  n_exceed <- sum(exceed_mask, na.rm = TRUE)
  prop_exceed <- if (total_cells_time > 0) n_exceed / total_cells_time else NA_real_
  # Conteo por tiempo (nº de celdas que exceden en cada t)
  exceed_by_time <- apply(exceed_mask, 3, function(m) sum(m, na.rm = TRUE))
  # Ratio BA/Área y peor caso
  ratio <- ba_arr / area_3d
  ratio[!is.finite(ratio)] <- NA
  max_ratio <- suppressWarnings(max(ratio, na.rm = TRUE))
  worst_idx <- if (is.finite(max_ratio)) which(ratio == max_ratio, arr.ind = TRUE)[1, ] else c(NA, NA, NA)
  
  worst <- if (all(is.finite(worst_idx))) {
    list(
      lon_idx = worst_idx[1],
      lat_idx = worst_idx[2],
      time_idx = worst_idx[3],
      lon = lon[worst_idx[1]],
      lat = lat[worst_idx[2]],
      BA = ba_arr[worst_idx[1], worst_idx[2], worst_idx[3]],
      Area = area_3d[worst_idx[1], worst_idx[2], worst_idx[3]],
      ratio = max_ratio
    )
  } else NULL
  
  # Top-10 excedencias con detalles
  w <- which(exceed_mask, arr.ind = TRUE)
  top10 <- NULL
  if (!is.null(w) && nrow(w) > 0) {
    ord <- order(ratio[w], decreasing = TRUE)
    top_n <- head(ord, 10)
    w_top <- w[top_n, , drop = FALSE]
    top10 <- data.frame(
      lon_idx = w_top[, 1],
      lat_idx = w_top[, 2],
      time_idx = w_top[, 3],
      lon = lon[w_top[, 1]],
      lat = lat[w_top[, 2]],
      BA = ba_arr[w_top],
      Area = area_3d[w_top],
      ratio = ratio[w_top]
    )
  }
  
  # Mapa de cuántas veces excede cada celda a lo largo del tiempo
  exceed_count_map <- apply(exceed_mask, c(1, 2), function(m) sum(m, na.rm = TRUE))
  list(
    any_exceed = any_exceed,
    n_exceed = n_exceed,
    total_cells_time = total_cells_time,
    prop_exceed = prop_exceed,
    exceed_by_time = exceed_by_time,
    worst = worst,
    top10 = top10,
    exceed_count_map = exceed_count_map,
    exceed_mask = exceed_mask # por si quieres inspeccionarlo/plotear
  )
}

# --- EJECUCIÓN ---
# Usa una tolerancia mínima para evitar falsos positivos por redondeo numérico.
res <- check_ba_exceeds_area(ba_arr = BA_qm_global_co,
                             area_mat = area_matrix, lon = lon, lat = lat, tol = 1e-6)
# --- RESÚMENES ÚTILES ---
cat("¿Hay excedencias? ", res$any_exceed, "\n")
cat("Excedencias totales (celda-tiempo): ", res$n_exceed, " de ", res$total_cells_time,
    sprintf(" (%.4f%%)\n", 100 * res$prop_exceed))
if (!is.null(res$worst)) {
  with(res$worst, {
    cat("Peor caso -> lon=", lon, ", lat=", lat, ", t_idx=", time_idx,
        ", BA=", BA, " km², Área=", Area, " km², Ratio=", ratio, "\n", sep = "")
  })
}
print(res$exceed_by_time)
if (!is.null(res$top10)) print(res$top10)



# ---------------------------
# TOPE FÍSICO MÁXIMO (lat-dependiente): BA ≤ area_matrix
# ---------------------------
ncols  <- dim(BA_qm_global_co)[1]
nrows  <- dim(BA_qm_global_co)[2]
ntime  <- dim(BA_qm_global_co)[3]

eps <- 1e-12  # tolerancia numérica
mask_sobre_max <- array(FALSE, dim = c(ncols, nrows, ntime))
for (k in seq_len(ntime)) {
  slice_k <- BA_qm_global_co[,,k]
  m <- (slice_k > (area_matrix + eps)) & is.finite(slice_k) & is.finite(area_matrix)
  # Recortar al máximo físico (no poner NA ni 0, sino el tope)
  slice_k[m] <- area_matrix[m]
  BA_qm_global_co[,,k] <- slice_k
  mask_sobre_max[,,k] <- m
}
cat("Celdas-tiempo recortadas por máximo: ", sum(mask_sobre_max), "\n")








ncols <- length(lon)   # debería ser 1440
nrows <- length(lat)   # debería ser 720
amin_px_eq  <- 0.09
amin_by_row <- amin_px_eq * cos(lat * pi/180)        # length = nrows
amin_matrix <- t(matrix(amin_by_row, nrow = nrows, ncol = ncols, byrow = FALSE))
# image.plot(lon, lat, amin_matrix)
dev.off()


mask_bajo_umbral <- array(FALSE, dim = c(ncols, nrows, ntime))  # máscara 3D
for (k in seq_len(ntime)) {
  slice_k <- BA_qm_global_co[,,k]
  m <- (slice_k < amin_matrix) & !is.na(slice_k) & !is.na(amin_matrix)  # máscara de cambios: valor real < umbral y no-NA
  slice_k[m] <- 0  # aplicar cambios
  BA_qm_global_co[,,k] <- slice_k  # escribir de vuelta
  mask_bajo_umbral[,,k] <- m
}
BA_qm_global_co[mask_bajo_umbral] <- 0
# image.plot(lon, lat, BA_harmonised[,,1])

rm(mask_bajo_umbral, mask_sobre_max, Mask_Fire51_aboveTope_global, 
   Mask_union_global, BA_Fire51_aboveTope_global, BA_FIRE60, slice_k); gc()
rm(amin_matrix, area_matrix, m, res, k, eps, cell_area_constant, nrows, nrow_grid, ntime, ncols, ncol_grid);gc()
rm(amin_by_row, amin_px_eq, area_by_row, t, check_ba_exceeds_area);gc()
BA_60=BA_qm_global_co
rm(BA_qm_global_co); gc()
BA_60_192<-BA_60[,,1:192]
BA_60_192[is.na(BA_60_192)]=0
rm(BA_60);gc()
BA_60_192=BA_60_192*1e6
# BA_60_192=BA_60[,,1:192]
rm(FIRECCI60_FINAL);gc()
# load(paste0(output_dir_RData, "FireCCIS311_S3_BA_monthly_2019_2024.RData"))
load(file.path(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))

dim(BA_60_192)
dim(s3)


stopifnot(identical(dim(BA_60_192)[1:2], dim(s3)[1:2]))

nxy <- dim(BA_60_192)[1:2]
out <- array(NA_real_, dim = c(nxy, 192 + 72))  # 1440 x 720 x 240
dim(out)
out[,, 1:192]   <- BA_60_192
out[,, 193:264] <- s3
# out[is.na(out)]=0
# Chequeo
dim(out)        # 1440 720 240
# sum(is.na(out)) # opcional
# (Opcional) vista rápida para verificar orientación
# image.plot(lon, lat, out[,,264], main = "BA_S3 (mes 1) - lat reorientada")

# output_dir_RData="/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-3/RData/"
BA_MRBA60=out
# --- Guardar a disco ---
out_rdata <- file.path(output_dir_RData, "MRBA60_BA_m2_monthly_2003_2024.RData")
save(BA_MRBA60, lon, lat, file = out_rdata)
cat("💾 Guardado:", out_rdata, "\n")
dim(BA_MRBA60)
