# =============================================================================
# EQM POSITIVOS-ONLY POR BIOMA x MES PARA MRBA60
# =============================================================================
# Objetivo:
# 1. Cargar FireCCI51, FireCCIS311 y BA armonizado.
# 2. Crear máscara de celdas con información.
# 3. Aplicar límites físicos mínimos y máximos.
# 4. Ajustar Quantile Mapping positivos-only por bioma x mes.
# 5. Aplicar EQM a 2003–2024.
# 6. Reponer valores FireCCI51 por encima del tope en 2003–2018.
# 7. Generar auditoría, Q-Q plots, histogramas y scatter.
# =============================================================================

rm(list = ls())
graphics.off()
gc()

# =============================================================================
# 1) LIBRERÍAS
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(lubridate)
  library(terra)
  library(raster)
  library(ggplot2)
  library(ncdf4)
  library(sp)
  library(fields)
  library(rworldmap)
  library(graticule)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(viridis)
  library(caret)
  library(randomForest)
  library(tidyr)
  library(tibble)
  library(cowplot)
  library(fastshap)
  library(qmap)
  library(scales)
  library(readr)
})

# =============================================================================
# 2) CONFIGURACIÓN GENERAL
# =============================================================================

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir          <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv      <- file.path(output_dir, "csv")
output_dir_plot     <- file.path(output_dir, "plot")
output_dir_RData    <- file.path(output_dir, "RData")
output_dir_plot_qm  <- file.path(output_dir, "plot_QM001M")
output_dir_qq       <- file.path(output_dir_plot, "qqplots_eqm")

dirs <- c(
  output_dir,
  output_dir_csv,
  output_dir_plot,
  output_dir_RData,
  output_dir_plot_qm,
  output_dir_qq
)

for (d in dirs) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# 3) PARÁMETROS DEL AJUSTE
# =============================================================================

# 2003-01 = índice 1
# 2019-01 = índice 193
# 2024-12 = índice 264
slice_train <- 193:264      # 2019–2024
thr_min     <- 0            # 0 = sin umbral adicional tras QM

min_cells_qm  <- 30         # mínimo de celdas positivas por año
min_years_pos <- 2          # mínimo de años válidos para activar QM

amin_px_eq <- 0.09          # mínimo físico ecuatorial, km2

# =============================================================================
# 4) FUNCIONES AUXILIARES
# =============================================================================

safe_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("[^[:alnum:]_]", "", x)
  x
}

wasserstein1d <- function(a, b, probs = seq(0, 1, by = 0.01)) {
  a <- a[is.finite(a) & a > 0]
  b <- b[is.finite(b) & b > 0]
  
  if (length(a) < 30 || length(b) < 30) return(NA_real_)
  
  Qa <- as.numeric(quantile(a, probs = probs, na.rm = TRUE, type = 7))
  Qb <- as.numeric(quantile(b, probs = probs, na.rm = TRUE, type = 7))
  
  mean(abs(Qa - Qb))
}

win3 <- function(m) {
  ((m + c(-1, 0, 1) - 1) %% 12) + 1
}

apply_min_threshold <- function(ba_arr, amin_matrix) {
  
  ntime <- dim(ba_arr)[3]
  mask_bajo_umbral <- array(FALSE, dim = dim(ba_arr))
  
  for (k in seq_len(ntime)) {
    slice_k <- ba_arr[, , k]
    
    m <- is.finite(slice_k) &
      is.finite(amin_matrix) &
      slice_k > 0 &
      slice_k < amin_matrix
    
    slice_k[m] <- 0
    ba_arr[, , k] <- slice_k
    mask_bajo_umbral[, , k] <- m
  }
  
  list(
    ba = ba_arr,
    mask = mask_bajo_umbral
  )
}

apply_max_threshold <- function(ba_arr, area_matrix, eps = 1e-12) {
  
  ntime <- dim(ba_arr)[3]
  mask_sobre_max <- array(FALSE, dim = dim(ba_arr))
  
  for (k in seq_len(ntime)) {
    slice_k <- ba_arr[, , k]
    
    m <- is.finite(slice_k) &
      is.finite(area_matrix) &
      slice_k > area_matrix + eps
    
    slice_k[m] <- area_matrix[m]
    ba_arr[, , k] <- slice_k
    mask_sobre_max[, , k] <- m
  }
  
  list(
    ba = ba_arr,
    mask = mask_sobre_max
  )
}

check_ba_exceeds_area <- function(ba_arr, area_mat, lon, lat, tol = 1e-6) {
  
  stopifnot(length(dim(ba_arr)) == 3)
  
  d <- dim(ba_arr)
  stopifnot(identical(dim(area_mat), d[1:2]))
  stopifnot(length(lon) == d[1])
  stopifnot(length(lat) == d[2])
  
  area_3d <- array(rep(area_mat, d[3]), dim = d)
  
  exceed_mask <- is.finite(ba_arr) &
    is.finite(area_3d) &
    ba_arr > area_3d + tol
  
  total_cells_time <- sum(is.finite(ba_arr))
  n_exceed <- sum(exceed_mask, na.rm = TRUE)
  prop_exceed <- if (total_cells_time > 0) n_exceed / total_cells_time else NA_real_
  
  exceed_by_time <- apply(exceed_mask, 3, function(m) sum(m, na.rm = TRUE))
  
  ratio <- ba_arr / area_3d
  ratio[!is.finite(ratio)] <- NA
  
  max_ratio <- suppressWarnings(max(ratio, na.rm = TRUE))
  
  worst_idx <- if (is.finite(max_ratio)) {
    which(ratio == max_ratio, arr.ind = TRUE)[1, ]
  } else {
    c(NA, NA, NA)
  }
  
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
  } else {
    NULL
  }
  
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
  
  exceed_count_map <- apply(exceed_mask, c(1, 2), function(m) sum(m, na.rm = TRUE))
  
  list(
    any_exceed = any(exceed_mask, na.rm = TRUE),
    n_exceed = n_exceed,
    total_cells_time = total_cells_time,
    prop_exceed = prop_exceed,
    exceed_by_time = exceed_by_time,
    worst = worst,
    top10 = top10,
    exceed_count_map = exceed_count_map,
    exceed_mask = exceed_mask
  )
}

df_hist_from_arrays <- function(list_arrays_named,
                                filter_zeros = TRUE,
                                sample_n = NULL,
                                force_zero = TRUE) {
  
  dfs <- lapply(names(list_arrays_named), function(nm) {
    
    v <- as.vector(list_arrays_named[[nm]])
    
    if (force_zero) v[!is.finite(v)] <- 0
    if (filter_zeros) v <- v[v > 0]
    
    tibble(
      dataset = nm,
      value = v
    )
  })
  
  df <- bind_rows(dfs)
  
  if (!is.null(sample_n) && nrow(df) > sample_n) {
    set.seed(123)
    df <- dplyr::sample_n(df, sample_n)
  }
  
  df
}

plot_hist_panel <- function(df, title, subtitle) {
  
  ggplot(df, aes(x = value, colour = dataset)) +
    geom_density(linewidth = 1, adjust = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    scale_colour_manual(
      values = c(
        "FireCCIS311"   = "blue",
        "FireCCI51"     = "orange",
        "Harmonised"    = "red",
        "Harmonised_QM" = "brown4"
      )
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Burned area (km²)",
      y = "Percentage",
      colour = "Dataset"
    ) +
    theme_minimal(base_size = 12)
}

# =============================================================================
# 5) CARGA DE LON/LAT Y CALENDARIO
# =============================================================================

load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))

lon_range <- as.vector(lon)
lat_range <- as.vector(lat)

dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

dates_all <- dates_full
mons_all  <- month(dates_all)
years_all <- year(dates_all)

# Mallas en formato lat x lon, útiles para sf
lon_mat <- matrix(
  rep(lon_range, each = length(lat_range)),
  nrow = length(lat_range),
  ncol = length(lon_range),
  byrow = FALSE
)

lat_mat <- matrix(
  rep(lat_range, times = length(lon_range)),
  nrow = length(lat_range),
  ncol = length(lon_range),
  byrow = FALSE
)

# =============================================================================
# 6) CARGA DE PRODUCTOS DE ÁREA QUEMADA
# =============================================================================

cat("\nCargando FireCCIS311...\n")

load(file.path(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))

BA_FireS3 <- s3 / 1e6
BA_FireS3[BA_FireS3 == 0] <- NA
rm(s3)
gc()

cat("\nCargando FireCCI51...\n")

load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

BA_Fire51_tot <- f51 / 1e6
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
rm(f51)
gc()

nx <- dim(BA_Fire51_tot)[1]
ny <- dim(BA_Fire51_tot)[2]
nt <- dim(BA_Fire51_tot)[3]

stopifnot(nt == length(dates_all))

# FireCCIS311 en cubo completo 2003–2024
BA_FireS3_tot <- array(NA_real_, dim = c(nx, ny, nt))
BA_FireS3_tot[, , ind_common] <- BA_FireS3

rm(BA_FireS3)
gc()

# =============================================================================
# 7) CARGA DEL PRODUCTO ARMONIZADO
# =============================================================================

cat("\nCargando BA armonizado corregido por tope FireCCI51...\n")

load(file.path(output_dir_RData, "BA_harmonised_correctedByF51Tope_B1-MRBA60-2003-2024.RData"))

stopifnot(exists("BA_harmonised"))
stopifnot(all(dim(BA_harmonised) == c(nx, ny, nt)))

BA_harmonised[BA_harmonised == 0] <- NA

# =============================================================================
# 8) MÁSCARA DE CELDAS CON INFORMACIÓN
# =============================================================================
# Equivale a tu bucle triple, pero vectorizado.
# Una celda-tiempo entra en la máscara si alguno de los tres productos tiene BA > 0.

cat("\nConstruyendo máscara común de información...\n")

status_matrix2_tot <- (
  (!is.na(BA_harmonised)  & BA_harmonised  > 0) |
    (!is.na(BA_Fire51_tot) & BA_Fire51_tot > 0) |
    (!is.na(BA_FireS3_tot) & BA_FireS3_tot > 0)
)

status_matrix2_tot <- ifelse(status_matrix2_tot, 1, NA)

BA_Fire51_tot[status_matrix2_tot == 1 & is.na(BA_Fire51_tot)] <- 0
BA_FireS3_tot[status_matrix2_tot == 1 & is.na(BA_FireS3_tot)] <- 0
BA_harmonised[status_matrix2_tot == 1 & is.na(BA_harmonised)] <- 0

gc()

# =============================================================================
# 9) MÁSCARAS FÍSICAS MÍNIMA Y MÁXIMA
# =============================================================================

cat("\nCalculando límites físicos por celda...\n")

cell_area_constant <- (110.57 * 0.25) * (111.32 * 0.25)

area_by_row <- cell_area_constant * cos(lat * pi / 180)

area_matrix <- matrix(
  area_by_row,
  nrow = length(lat),
  ncol = length(lon),
  byrow = FALSE
)

area_matrix <- t(area_matrix)

stopifnot(identical(dim(area_matrix), c(nx, ny)))

amin_by_row <- amin_px_eq * cos(lat * pi / 180)

amin_matrix <- matrix(
  amin_by_row,
  nrow = length(lat),
  ncol = length(lon),
  byrow = FALSE
)

amin_matrix <- t(amin_matrix)

stopifnot(identical(dim(amin_matrix), c(nx, ny)))

# Aplicar mínimo físico al armonizado de entrada
tmp_min <- apply_min_threshold(BA_harmonised, amin_matrix)
BA_harmonised <- tmp_min$ba
mask_bajo_umbral_input <- tmp_min$mask

cat("Celdas-tiempo bajo umbral mínimo en BA_harmonised:",
    sum(mask_bajo_umbral_input), "\n")

rm(tmp_min)
gc()

# =============================================================================
# 10) CARGA DE BIOMAS
# =============================================================================

cat("\nCargando biomas...\n")

biomas_shp <- st_read(
  file.path(dir_oss, "continental-biomes_dinerstein_V10.shp"),
  quiet = TRUE
)

biomas_shp <- st_transform(biomas_shp, crs = 4326)

stopifnot("cont_bm" %in% names(biomas_shp))

# Asignación global rápida para obtener lista de biomas presentes
grid_points <- st_as_sf(
  data.frame(
    lon = as.vector(lon_mat),
    lat = as.vector(lat_mat)
  ),
  coords = c("lon", "lat"),
  crs = 4326
)

inter <- st_intersects(grid_points, biomas_shp)

bioma_asignado <- sapply(inter, function(i) {
  if (length(i) == 0) NA_character_ else biomas_shp$cont_bm[i[1]]
})

grid_points$bioma_final <- bioma_asignado

biomas_unique <- sort(unique(na.omit(grid_points$bioma_final)))

cat("Número de biomas detectados:", length(biomas_unique), "\n")

rm(grid_points, inter, bioma_asignado)
gc()

# =============================================================================
# 11) INICIALIZACIÓN DE SALIDA
# =============================================================================

BA_qm_global <- array(NA_real_, dim = dim(BA_harmonised))
audit_rows <- list()

# =============================================================================
# 12) LOOP PRINCIPAL POR BIOMA
# =============================================================================

for (bioma in biomas_unique) {
  
  cat("\nProcesando bioma:", bioma, "\n")
  
  safe_biome_name <- safe_name(bioma)
  
  bioma_sel <- biomas_shp %>% filter(cont_bm == bioma)
  bbox <- st_bbox(bioma_sel)
  
  lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
  lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])
  
  if (length(lon_idx) == 0 || length(lat_idx) == 0) {
    cat("  Sin celdas en bbox. Se omite.\n")
    next
  }
  
  lon_mat_crop <- lon_mat[lat_idx, lon_idx, drop = FALSE]
  lat_mat_crop <- lat_mat[lat_idx, lon_idx, drop = FALSE]
  
  status_crop <- status_matrix2_tot[lon_idx, lat_idx, , drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # 12.1) Máscara del bioma en el recorte
  # ---------------------------------------------------------------------------
  
  grid_points_crop <- st_as_sf(
    data.frame(
      lon = as.vector(lon_mat_crop),
      lat = as.vector(lat_mat_crop)
    ),
    coords = c("lon", "lat"),
    crs = 4326
  )
  
  inter_crop <- st_intersects(grid_points_crop, biomas_shp)
  
  bioma_crop <- sapply(inter_crop, function(i) {
    if (length(i) == 0) "Ninguno" else biomas_shp$cont_bm[i[1]]
  })
  
  grid_points_crop$bioma_final <- bioma_crop
  
  presencia_en_serie <- apply(status_crop, c(1, 2), function(x) {
    any(x == 1, na.rm = TRUE)
  })
  
  grid_points_crop$presencia_serie <- as.vector(t(presencia_en_serie))
  
  puntos_sin_bioma_con_presencia <- grid_points_crop %>%
    filter(bioma_final == "Ninguno" & presencia_serie)
  
  puntos_con_bioma <- grid_points_crop %>%
    filter(bioma_final != "Ninguno")
  
  if (nrow(puntos_sin_bioma_con_presencia) > 0 && nrow(puntos_con_bioma) > 0) {
    
    nearest_idx <- st_nearest_feature(
      puntos_sin_bioma_con_presencia,
      puntos_con_bioma
    )
    
    idx_replace <- which(
      grid_points_crop$bioma_final == "Ninguno" &
        grid_points_crop$presencia_serie
    )
    
    grid_points_crop$bioma_final[idx_replace] <-
      puntos_con_bioma$bioma_final[nearest_idx]
  }
  
  mask_vec <- grid_points_crop$bioma_final == bioma
  
  mask_latlon <- matrix(
    mask_vec,
    nrow = length(lat_idx),
    ncol = length(lon_idx),
    byrow = FALSE
  )
  
  mask_final_clean <- t(mask_latlon)
  
  stopifnot(identical(dim(mask_final_clean), c(length(lon_idx), length(lat_idx))))
  
  if (!any(mask_final_clean, na.rm = TRUE)) {
    cat("  Máscara vacía. Se omite.\n")
    next
  }
  
  # ---------------------------------------------------------------------------
  # 12.2) Recortes y enmascarado del bioma
  # ---------------------------------------------------------------------------
  
  BA_FireS3_crop_tot <- BA_FireS3_tot[lon_idx, lat_idx, , drop = FALSE]
  BA_harm_crop_tot   <- BA_harmonised[lon_idx, lat_idx, , drop = FALSE]
  
  for (tt in seq_len(nt)) {
    
    s3 <- BA_FireS3_crop_tot[, , tt]
    hm <- BA_harm_crop_tot[, , tt]
    
    s3[!mask_final_clean] <- NA
    hm[!mask_final_clean] <- NA
    
    BA_FireS3_crop_tot[, , tt] <- s3
    BA_harm_crop_tot[, , tt]   <- hm
  }
  
  # ---------------------------------------------------------------------------
  # 12.3) Si no hay S3 positivo en 2019–2024, copiar armonizado con límites
  # ---------------------------------------------------------------------------
  
  if (sum(BA_FireS3_crop_tot[, , slice_train], na.rm = TRUE) == 0) {
    
    cat("  Sin BA positiva en FireCCIS311 durante calibración. Se copia armonizado.\n")
    
    max_mask_crop <- area_matrix[lon_idx, lat_idx, drop = FALSE]
    min_mask_crop <- amin_matrix[lon_idx, lat_idx, drop = FALSE]
    
    for (tt in seq_len(nt)) {
      
      slice_out <- BA_harm_crop_tot[, , tt]
      
      slice_out[slice_out > 0 & slice_out <= min_mask_crop] <- 0
      slice_out <- pmin(slice_out, max_mask_crop)
      
      BA_harm_crop_tot[, , tt] <- slice_out
    }
    
    for (tt in seq_len(nt)) {
      tmp_global <- BA_qm_global[lon_idx, lat_idx, tt]
      tmp_global[mask_final_clean] <- BA_harm_crop_tot[, , tt][mask_final_clean]
      BA_qm_global[lon_idx, lat_idx, tt] <- tmp_global
    }
    
    rm(BA_FireS3_crop_tot, BA_harm_crop_tot)
    gc()
    
    next
  }
  
  # ---------------------------------------------------------------------------
  # 12.4) Entrenamiento QM por mes, usando ventana móvil de 3 meses
  # ---------------------------------------------------------------------------
  
  train_years <- sort(unique(years_all[slice_train]))
  train_cache <- vector("list", 12)
  
  for (m in 1:12) {
    
    t_train_win <- integer(0)
    
    for (yy in train_years) {
      
      mset <- win3(m)
      
      ti <- which(years_all == yy & mons_all %in% mset)
      ti <- intersect(ti, slice_train)
      
      if (length(ti) > 0) {
        t_train_win <- c(t_train_win, ti)
      }
    }
    
    t_train_win <- sort(unique(t_train_win))
    
    if (length(t_train_win) == 0) {
      train_cache[[m]] <- list(apply_qm = FALSE, fit = NULL)
      next
    }
    
    ok_years <- 0
    
    for (yy in train_years) {
      
      ti_y <- intersect(t_train_win, which(years_all == yy))
      
      if (length(ti_y) == 0) next
      
      n_pos_yy <- sum(BA_FireS3_crop_tot[, , ti_y] > 0, na.rm = TRUE)
      
      if (n_pos_yy >= min_cells_qm) ok_years <- ok_years + 1
    }
    
    if (ok_years < min_years_pos) {
      train_cache[[m]] <- list(apply_qm = FALSE, fit = NULL)
      next
    }
    
    mask_rep <- rep(mask_final_clean, length(t_train_win))
    
    v_est <- as.vector(BA_harm_crop_tot[, , t_train_win])[mask_rep]
    v_ref <- as.vector(BA_FireS3_crop_tot[, , t_train_win])[mask_rep]
    
    keep <- is.finite(v_est) & is.finite(v_ref)
    
    v_est <- v_est[keep]
    v_ref <- v_ref[keep]
    
    est_pos <- v_est[v_est > 0]
    ref_pos <- v_ref[v_ref > 0]
    
    if (
      length(est_pos) < 30 ||
      length(ref_pos) < 30 ||
      length(unique(est_pos)) < 30 ||
      length(unique(ref_pos)) < 30
    ) {
      train_cache[[m]] <- list(apply_qm = FALSE, fit = NULL)
      next
    }
    
    fit <- try(
      qmap::fitQmapQUANT(
        obs = ref_pos,
        mod = est_pos,
        qstep = 0.01,
        type = "linear"
      ),
      silent = TRUE
    )
    
    if (inherits(fit, "try-error")) {
      
      train_cache[[m]] <- list(apply_qm = FALSE, fit = NULL)
      
    } else {
      
      train_cache[[m]] <- list(apply_qm = TRUE, fit = fit)
      
      cat(sprintf(
        "  [TRAIN-3M] %s | mes=%02d | n_ref_pos=%d | n_est_pos=%d | años_ok=%d/%d\n",
        safe_biome_name,
        m,
        length(ref_pos),
        length(est_pos),
        ok_years,
        length(train_years)
      ))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 12.5) Aplicación QM a todo 2003–2024
  # ---------------------------------------------------------------------------
  
  BA_qm_crop_tot <- array(NA_real_, dim = dim(BA_harm_crop_tot))
  
  max_mask_crop <- area_matrix[lon_idx, lat_idx, drop = FALSE]
  min_mask_crop <- amin_matrix[lon_idx, lat_idx, drop = FALSE]
  
  for (m in 1:12) {
    
    t_apply <- which(mons_all == m)
    tr <- train_cache[[m]]
    
    for (tt in t_apply) {
      
      slice_est <- BA_harm_crop_tot[, , tt]
      slice_out <- slice_est
      
      vals <- slice_est[mask_final_clean]
      keep_apply <- is.finite(vals)
      x <- vals[keep_apply]
      
      if (!is.null(tr) && isTRUE(tr$apply_qm) && !is.null(tr$fit)) {
        
        pos_idx <- which(is.finite(x) & x > 0)
        
        if (length(pos_idx) > 0) {
          
          mapped <- try(
            qmap::doQmapQUANT(
              x[pos_idx],
              tr$fit,
              type = "linear"
            ),
            silent = TRUE
          )
          
          if (!inherits(mapped, "try-error")) {
            
            if (thr_min > 0) {
              mapped[mapped > 0 & mapped < thr_min] <- 0
            }
            
            mapped[mapped < 0] <- 0
            x[pos_idx] <- mapped
          }
        }
      }
      
      vals[keep_apply] <- x
      slice_out[mask_final_clean] <- vals
      
      slice_out[slice_out > 0 & slice_out <= min_mask_crop] <- 0
      slice_out <- pmin(slice_out, max_mask_crop)
      
      BA_qm_crop_tot[, , tt] <- slice_out
    }
  }
  
  # ---------------------------------------------------------------------------
  # 12.6) Auditoría 2019–2024
  # ---------------------------------------------------------------------------
  
  for (m in 1:12) {
    
    t_cal <- slice_train[mons_all[slice_train] == m]
    
    if (length(t_cal) == 0) next
    
    est_fold <- c()
    ref_fold <- c()
    adj_fold <- c()
    
    for (tt in t_cal) {
      
      ve0 <- BA_harm_crop_tot[, , tt][mask_final_clean]
      vr  <- BA_FireS3_crop_tot[, , tt][mask_final_clean]
      ve1 <- BA_qm_crop_tot[, , tt][mask_final_clean]
      
      est_fold <- c(est_fold, ve0[is.finite(ve0) & ve0 > 0])
      ref_fold <- c(ref_fold, vr[is.finite(vr) & vr > 0])
      adj_fold <- c(adj_fold, ve1[is.finite(ve1) & ve1 > 0])
    }
    
    rmse_before <- NA_real_
    rmse_after  <- NA_real_
    mae_before  <- NA_real_
    mae_after   <- NA_real_
    
    n0 <- min(length(est_fold), length(ref_fold))
    n1 <- min(length(adj_fold), length(ref_fold))
    
    if (n0 >= 30) {
      set.seed(123)
      ei <- sample(seq_along(est_fold), n0)
      ri <- sample(seq_along(ref_fold), n0)
      
      rmse_before <- sqrt(mean((est_fold[ei] - ref_fold[ri])^2))
      mae_before  <- mean(abs(est_fold[ei] - ref_fold[ri]))
    }
    
    if (n1 >= 30) {
      set.seed(123)
      ai <- sample(seq_along(adj_fold), n1)
      ri <- sample(seq_along(ref_fold), n1)
      
      rmse_after <- sqrt(mean((adj_fold[ai] - ref_fold[ri])^2))
      mae_after  <- mean(abs(adj_fold[ai] - ref_fold[ri]))
    }
    
    bias_mean_before <- if (length(est_fold) >= 30 && length(ref_fold) >= 30) {
      mean(est_fold, na.rm = TRUE) - mean(ref_fold, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    bias_mean_after <- if (length(adj_fold) >= 30 && length(ref_fold) >= 30) {
      mean(adj_fold, na.rm = TRUE) - mean(ref_fold, na.rm = TRUE)
    } else {
      NA_real_
    }
    
    probs <- c(0.10, 0.50, 0.90)
    
    q_est <- if (length(est_fold) >= 30) {
      as.numeric(quantile(est_fold, probs, na.rm = TRUE))
    } else {
      rep(NA_real_, 3)
    }
    
    q_adj <- if (length(adj_fold) >= 30) {
      as.numeric(quantile(adj_fold, probs, na.rm = TRUE))
    } else {
      rep(NA_real_, 3)
    }
    
    q_ref <- if (length(ref_fold) >= 30) {
      as.numeric(quantile(ref_fold, probs, na.rm = TRUE))
    } else {
      rep(NA_real_, 3)
    }
    
    w1_before <- wasserstein1d(est_fold, ref_fold)
    w1_after  <- wasserstein1d(adj_fold, ref_fold)
    
    dW1 <- if (is.finite(w1_before) && is.finite(w1_after)) {
      w1_before - w1_after
    } else {
      NA_real_
    }
    
    # Q-Q plot
    if (length(ref_fold) >= 30 && (length(est_fold) >= 30 || length(adj_fold) >= 30)) {
      
      qgrid <- seq(0, 1, by = 0.01)
      
      q_ref <- as.numeric(quantile(ref_fold, qgrid, na.rm = TRUE))
      
      q_est <- if (length(est_fold) >= 30) {
        as.numeric(quantile(est_fold, qgrid, na.rm = TRUE))
      } else {
        rep(NA_real_, length(qgrid))
      }
      
      q_adj <- if (length(adj_fold) >= 30) {
        as.numeric(quantile(adj_fold, qgrid, na.rm = TRUE))
      } else {
        rep(NA_real_, length(qgrid))
      }
      
      dfqq <- data.frame(
        q = qgrid,
        ref = q_ref,
        est = q_est,
        adj = q_adj
      )
      
      pqq <- ggplot(dfqq, aes(x = ref)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        { if (any(is.finite(dfqq$est))) geom_path(aes(y = est), alpha = 0.6) } +
        { if (any(is.finite(dfqq$adj))) geom_path(aes(y = adj), alpha = 0.8, color = "steelblue") } +
        labs(
          title = paste0("Q-Q vs S3 | ", safe_biome_name, " | month ", sprintf("%02d", m)),
          x = "S3 quantiles",
          y = "Quantiles: before = black, after = blue"
        ) +
        theme_minimal(base_size = 11)
      
      ggsave(
        file.path(
          output_dir_qq,
          paste0("qq_", safe_biome_name, "_m", sprintf("%02d", m), ".png")
        ),
        pqq,
        width = 6.5,
        height = 5.2,
        dpi = 220
      )
    }
    
    audit_rows[[length(audit_rows) + 1]] <- tibble(
      biome = safe_biome_name,
      month = m,
      n_pos_ref_calib = length(ref_fold),
      n_pos_est_calib = length(est_fold),
      n_pos_adj_calib = length(adj_fold),
      RMSE_before = rmse_before,
      RMSE_after = rmse_after,
      dRMSE = ifelse(
        is.na(rmse_before) | is.na(rmse_after),
        NA_real_,
        rmse_before - rmse_after
      ),
      MAE_before = mae_before,
      MAE_after = mae_after,
      dMAE = ifelse(
        is.na(mae_before) | is.na(mae_after),
        NA_real_,
        mae_before - mae_after
      ),
      BiasMean_before = bias_mean_before,
      BiasMean_after = bias_mean_after,
      BiasP10_before = q_est[1] - q_ref[1],
      BiasP10_after = q_adj[1] - q_ref[1],
      BiasP50_before = q_est[2] - q_ref[2],
      BiasP50_after = q_adj[2] - q_ref[2],
      BiasP90_before = q_est[3] - q_ref[3],
      BiasP90_after = q_adj[3] - q_ref[3],
      W1_before = w1_before,
      W1_after = w1_after,
      dW1 = dW1
    )
  }
  
  # ---------------------------------------------------------------------------
  # 12.7) Volcado al cubo global
  # ---------------------------------------------------------------------------
  
  for (tt in seq_len(nt)) {
    
    biome_values <- BA_qm_crop_tot[, , tt]
    tmp_global <- BA_qm_global[lon_idx, lat_idx, tt]
    
    tmp_global[mask_final_clean] <- biome_values[mask_final_clean]
    
    BA_qm_global[lon_idx, lat_idx, tt] <- tmp_global
  }
  
  rm(
    grid_points_crop,
    inter_crop,
    bioma_crop,
    presencia_en_serie,
    BA_FireS3_crop_tot,
    BA_harm_crop_tot,
    BA_qm_crop_tot
  )
  
  gc()
}

# =============================================================================
# 13) GUARDAR AUDITORÍA
# =============================================================================
# 
# audit_df <- bind_rows(audit_rows)
# 
# audit_file <- file.path(
#   output_dir_csv,
#   paste0("AUDIT_EQM_positive_only_", Modelo, "_2019_2024.csv")
# )
# 
# write_csv(audit_df, audit_file)
# 
# cat("\nAuditoría guardada en:\n", audit_file, "\n")

# =============================================================================
# 14) HISTOGRAMAS / DENSIDADES 2019–2024
# =============================================================================

sel_range <- 193:264

BA_FireS3_sel  <- BA_FireS3_tot[, , sel_range]
BA_Fire51_sel  <- BA_Fire51_tot[, , sel_range]
BA_harm_sel    <- BA_harmonised[, , sel_range]
BA_qm_sel      <- BA_qm_global[, , sel_range]

# Panel A: FireCCI51 vs FireCCIS311
hist_A <- df_hist_from_arrays(list(
  FireCCIS311 = BA_FireS3_sel,
  FireCCI51   = BA_Fire51_sel
))

p_hist_A <- plot_hist_panel(
  hist_A,
  "Panel a): FireCCI51 & FireCCIS311",
  "Distribution, 2019–2024"
)

ggsave(
  file.path(output_dir_plot_qm, "Panel_a_Hist_FireS3_vs_FireCCI51.png"),
  p_hist_A,
  width = 10,
  height = 8,
  dpi = 300
)

# Panel B: añade armonizado antes de QM
hist_B <- df_hist_from_arrays(list(
  FireCCIS311 = BA_FireS3_sel,
  FireCCI51   = BA_Fire51_sel,
  Harmonised  = BA_harm_sel
))

p_hist_B <- plot_hist_panel(
  hist_B,
  "Panel b): FireCCI51, FireCCIS311 & Harmonised",
  "Distribution, 2019–2024"
)

ggsave(
  file.path(output_dir_plot_qm, "Panel_b_Hist_FireS3_FireCCI51_Harmonised.png"),
  p_hist_B,
  width = 10,
  height = 8,
  dpi = 300
)

# Panel C: añade armonizado tras QM
hist_C <- df_hist_from_arrays(list(
  FireCCIS311   = BA_FireS3_sel,
  FireCCI51     = BA_Fire51_sel,
  Harmonised    = BA_harm_sel,
  Harmonised_QM = BA_qm_sel
))

p_hist_C <- plot_hist_panel(
  hist_C,
  "Panel c): FireCCI51, FireCCIS311, Harmonised & Harmonised QM",
  "Distribution, 2019–2024"
)

ggsave(
  file.path(output_dir_plot_qm, "Panel_c_Hist_FireS3_FireCCI51_Harmonised_QM.png"),
  p_hist_C,
  width = 10,
  height = 8,
  dpi = 300
)

rm(
  BA_FireS3_sel,
  BA_Fire51_sel,
  BA_harm_sel,
  BA_qm_sel,
  hist_A,
  hist_B,
  hist_C,
  p_hist_A,
  p_hist_B,
  p_hist_C
)

gc()

# =============================================================================
# 15) REPOSICIÓN DE VALORES FIRECCI51 POR ENCIMA DEL TOPE EN 2003–2018
# =============================================================================

cat("\nCargando máscaras de FireCCI51 por encima del tope...\n")

load(file.path(
  output_dir_RData,
  paste0("F51_maskAboveTope_GLOBAL_", Modelo, "_2003_2024.RData")
))

load(file.path(
  output_dir_RData,
  paste0("F51_aboveTope_GLOBAL_", Modelo, "_2003_2024.RData")
))

Mask_Fire51_aboveTope_global <- 
  (Mask_Fire51_aboveTope_global == 1) |
  (Mask_Fire51_aboveTope_global == TRUE)

Mask_Fire51_aboveTope_global[is.na(Mask_Fire51_aboveTope_global)] <- FALSE

BA_qm_global_co <- BA_qm_global

for (tt in 1:192) {
  
  m <- Mask_Fire51_aboveTope_global[, , tt]
  
  if (!any(m, na.rm = TRUE)) next
  
  BA_qm_global_co[, , tt][m] <- BA_Fire51_aboveTope_global[, , tt][m]
}

rm(Mask_Fire51_aboveTope_global, BA_Fire51_aboveTope_global)
gc()

# =============================================================================
# 16) COMPROBACIÓN Y APLICACIÓN DE LÍMITE FÍSICO MÁXIMO
# =============================================================================

cat("\nComprobando excedencias BA > área de celda antes del recorte...\n")

res_exceed <- check_ba_exceeds_area(
  ba_arr = BA_qm_global_co,
  area_mat = area_matrix,
  lon = lon,
  lat = lat,
  tol = 1e-6
)

cat("¿Hay excedencias?:", res_exceed$any_exceed, "\n")
cat(
  "Excedencias totales:",
  res_exceed$n_exceed,
  "de",
  res_exceed$total_cells_time,
  sprintf("(%.4f%%)\n", 100 * res_exceed$prop_exceed)
)

if (!is.null(res_exceed$worst)) {
  with(res_exceed$worst, {
    cat(
      "Peor caso -> lon=", lon,
      ", lat=", lat,
      ", t_idx=", time_idx,
      ", BA=", BA,
      " km², Área=", Area,
      " km², Ratio=", ratio,
      "\n",
      sep = ""
    )
  })
}

if (!is.null(res_exceed$top10)) {
  print(res_exceed$top10)
}

tmp_max <- apply_max_threshold(
  BA_qm_global_co,
  area_matrix,
  eps = 1e-12
)

BA_qm_global_co <- tmp_max$ba
mask_sobre_max <- tmp_max$mask

cat("Celdas-tiempo recortadas por máximo físico:",
    sum(mask_sobre_max), "\n")

rm(tmp_max)
gc()

# =============================================================================
# 17) APLICACIÓN FINAL DEL UMBRAL MÍNIMO
# =============================================================================

tmp_min_final <- apply_min_threshold(
  BA_qm_global_co,
  amin_matrix
)

BA_qm_global_co <- tmp_min_final$ba
mask_bajo_umbral_final <- tmp_min_final$mask

cat("Celdas-tiempo puestas a 0 por mínimo físico final:",
    sum(mask_bajo_umbral_final), "\n")

rm(tmp_min_final)
gc()


# =============================================================================
# 19) SUMAS ANUALES 2003–2018
# =============================================================================

cat("\nCalculando acumulados anuales 2003–2018...\n")

BA_51_sub <- BA_Fire51_tot[, , 1:192]
BA_harmonised_sub <- BA_qm_global_co[, , 1:192]

n_years <- 16
months_per_year <- 12

BA_51_annual <- array(
  NA_real_,
  dim = c(dim(BA_51_sub)[1], dim(BA_51_sub)[2], n_years)
)

BA_harmonised_annual <- array(
  NA_real_,
  dim = c(dim(BA_harmonised_sub)[1], dim(BA_harmonised_sub)[2], n_years)
)

for (yy in 1:n_years) {
  
  idx <- ((yy - 1) * months_per_year + 1):(yy * months_per_year)
  
  BA_51_annual[, , yy] <- apply(
    BA_51_sub[, , idx, drop = FALSE],
    c(1, 2),
    sum,
    na.rm = TRUE
  )
  
  BA_harmonised_annual[, , yy] <- apply(
    BA_harmonised_sub[, , idx, drop = FALSE],
    c(1, 2),
    sum,
    na.rm = TRUE
  )
}

rm(BA_51_sub, BA_harmonised_sub)
gc()

# =============================================================================
# 20) HISTOGRAMA FINAL CON PRODUCTO CORREGIDO
# =============================================================================

cat("\nGenerando histograma final con Harmonised_QM_corrected...\n")

BA_FireS3_sel <- BA_FireS3_tot[, , sel_range]
BA_Fire51_sel <- BA_Fire51_tot[, , sel_range]
BA_harm_sel   <- BA_harmonised[, , sel_range]
BA_qm_co_sel  <- BA_qm_global_co[, , sel_range]

hist_final <- df_hist_from_arrays(list(
  FireCCIS311   = BA_FireS3_sel,
  FireCCI51     = BA_Fire51_sel,
  Harmonised    = BA_harm_sel,
  Harmonised_QM = BA_qm_co_sel
))

p_hist_final <- plot_hist_panel(
  hist_final,
  "FireCCI51, FireCCIS311, Harmonised & Harmonised QM corrected",
  "Distribution, 2019–2024"
)

ggsave(
  file.path(output_dir_plot_qm, "Panel_c_Hist_FireS3_FireCCI51_Harmonised_QM_co.png"),
  p_hist_final,
  width = 10,
  height = 8,
  dpi = 300
)

rm(
  BA_FireS3_sel,
  BA_Fire51_sel,
  BA_harm_sel,
  BA_qm_co_sel,
  hist_final,
  p_hist_final
)

gc()



sel_range <- 1:192
BA_qm_sel <- BA_qm_global_co[,,sel_range]
BA_harm_sel <- BA_harmonised[,,sel_range]
BA_Fire51_sel <- BA_Fire51_tot[,,sel_range]
output_dir_plot_rle <- paste0(output_dir, "/plot_QM001M/")
library(scales)



# --- Helper: density curves with percentage on Y axis ---
plot_hist_panel <- function(df, title, subtitle) {
  ggplot(df, aes(x = value, colour = dataset)) +
    geom_density(linewidth = 1, adjust = 1) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),     # 10^1, 10^2, 10^3...
      labels = trans_format("log10", math_format(10^.x))    # pretty labels
    ) +
    scale_colour_manual(
      values = c("FireCCI51" = "orange", "Harmonised" = "red",
                 "Harmonised_QM" = "brown4")
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Burned area (km²)",
      y = "Percentage",
      colour = "Dataset"
    ) +
    theme_minimal(base_size = 12)
}

# --- Example usage ---
hist_A <- df_hist_from_arrays(list(
  FireCCI51   = BA_Fire51_sel,
  Harmonised = BA_harm_sel,
  Harmonised_QM = BA_qm_sel
))

p_hist_A <- plot_hist_panel(
  hist_A,
  "Panel c): FireCCI51, Harmonised & Harmonised QM",
  "Distribution (2003–2018)"
)

p_hist_A
# --- Guardar figuras ---
ggsave(file.path(output_dir_plot_rle, "Panel_d_2003-2018_Hist_FireCCI51_Harmonised_QM.png"), p_hist_A, width = 10, height = 8, dpi = 300)


# # =============================================================================
# # 21) GUARDAR RESULTADO FINAL
# # =============================================================================
BA_FIRE60=BA_qm_global_co
outfile <- file.path(output_dir_RData, "BA_MRBA60.RData")
save(BA_FIRE60, file = outfile)

