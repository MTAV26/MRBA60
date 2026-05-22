
process_biome_selected_predictors <- function(bioma) {
t_start <- Sys.time()
safe_biome_name <- safe_name(as.character(bioma))

log_file <- file.path(output_dir_log, paste0("LOG_SelectedPredictors_", safe_biome_name, ".txt"))
log_con <- file(log_file, open = "wt")

sink(log_con)
sink(log_con, type = "message")

on.exit({
  sink(type = "message")
  sink()
  close(log_con)
}, add = TRUE)

cat("============================================================\n")
cat("Procesando selección de predictores para bioma:", bioma, "\n")
cat("Inicio:", as.character(t_start), "\n")
cat("PID:", Sys.getpid(), "\n")
cat("============================================================\n")

biome_id <- match(bioma, biomas_unique)

# ---------------------------------------------------------------------------
# 7.1 Selección espacial del bioma
# ---------------------------------------------------------------------------

bioma_sel <- biomas_shp %>%
  dplyr::filter(cont_bm == bioma)

bbox <- sf::st_bbox(bioma_sel)

lon_idx <- which(lon_range >= bbox["xmin"] & lon_range <= bbox["xmax"])
lat_idx <- which(lat_range >= bbox["ymin"] & lat_range <= bbox["ymax"])

if (length(lon_idx) == 0 || length(lat_idx) == 0) {
  cat("Sin celdas para el bioma:", bioma, "\n")
  
  return(list(
    Biome = safe_biome_name,
    Status = "SKIPPED_NO_GRID_CELLS",
    n_months_selected = 0L,
    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  ))
}

crop3 <- function(x) {
  x[lon_idx, lat_idx, , drop = FALSE]
}

status_crop <- crop3(status_matrix2_tot)

lon_mat_crop <- lon_mat[lat_idx, lon_idx, drop = FALSE]
lat_mat_crop <- lat_mat[lat_idx, lon_idx, drop = FALSE]

lon_vec <- lon[lon_idx]
lat_vec <- lat[lat_idx]

# ---------------------------------------------------------------------------
# 7.2 Asignación del bioma a la grilla
# ---------------------------------------------------------------------------

grid_points_crop <- sf::st_as_sf(
  data.frame(
    lon = as.vector(lon_mat_crop),
    lat = as.vector(lat_mat_crop)
  ),
  coords = c("lon", "lat"),
  crs = 4326
)

inter <- sf::st_intersects(grid_points_crop, biomas_shp)

bioma_asignado <- sapply(inter, function(i) {
  if (length(i) == 0) return("Ninguno")
  biomas_shp$cont_bm[i[1]]
})

grid_points_biomas <- grid_points_crop
grid_points_biomas$bioma_final <- bioma_asignado

presencia_en_serie <- apply(
  status_crop,
  c(1, 2),
  function(x) any(x == 1, na.rm = TRUE)
)

grid_points_biomas$presencia_serie <- as.vector(t(presencia_en_serie))

puntos_sin_bioma_con_presencia <- grid_points_biomas %>%
  dplyr::filter(bioma_final == "Ninguno" & presencia_serie)

puntos_con_bioma <- grid_points_biomas %>%
  dplyr::filter(bioma_final != "Ninguno")

if (
  nrow(puntos_sin_bioma_con_presencia) > 0 &&
  nrow(puntos_con_bioma) > 0
) {
  
  nearest_idx <- sf::st_nearest_feature(
    puntos_sin_bioma_con_presencia,
    puntos_con_bioma
  )
  
  grid_points_biomas$bioma_final[
    which(
      grid_points_biomas$bioma_final == "Ninguno" &
        grid_points_biomas$presencia_serie
    )
  ] <- puntos_con_bioma$bioma_final[nearest_idx]
}

mask <- grid_points_biomas$bioma_final == bioma

mask_matrix <- matrix(
  mask,
  nrow = length(lat_vec),
  ncol = length(lon_vec),
  byrow = FALSE
)

mask_lonlat <- t(mask_matrix)

if (!any(mask_lonlat, na.rm = TRUE)) {
  cat("Máscara vacía para bioma:", bioma, "\n")
  
  return(list(
    Biome = safe_biome_name,
    Status = "SKIPPED_EMPTY_MASK",
    n_months_selected = 0L,
    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  ))
}

# ---------------------------------------------------------------------------
# 7.3 Recorte y máscara de predictores
# ---------------------------------------------------------------------------

cat("Recortando arrays...\n")

arr <- list(
  s3 = crop3(BA_FireS3_tot),
  f51 = crop3(BA_Fire51_tot),
  af = crop3(count_ActiveFire_tot),
  prec = crop3(PREC_tot),
  temp = crop3(TEMP_tot),
  frp_sum = crop3(FRPsum_tot),
  frp_median = crop3(FRPmedian_tot),
  fwi = crop3(FWI_tot),
  ndvi = crop3(NDVI_tot),
  wind = crop3(WIND_tot),
  cloud = crop3(cloud_tot),
  vpd = crop3(VPD_tot),
  soil = crop3(SOIL_tot)
)

arr <- lapply(arr, apply_mask_3d_2, mask_lonlat = mask_lonlat)

lat_2d <- t(lat_mat_crop)
lon_2d <- t(lon_mat_crop)

lat_2d[!mask_lonlat] <- NA
lon_2d[!mask_lonlat] <- NA

arr$lat <- array(
  rep(lat_2d, ntime_full),
  dim = c(dim(lat_2d), ntime_full)
)

arr$lon <- array(
  rep(lon_2d, ntime_full),
  dim = c(dim(lon_2d), ntime_full)
)

arr_common <- lapply(
  arr,
  function(x) x[, , ind_common, drop = FALSE]
)

total_fire <- sum(arr_common$s3, na.rm = TRUE)

if (!is.finite(total_fire) || total_fire == 0) {
  cat("El bioma no registra incendios en el periodo common. Se omite.\n")
  
  return(list(
    Biome = safe_biome_name,
    Status = "SKIPPED_NO_COMMON_FIRE",
    n_months_selected = 0L,
    seconds = as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  ))
}

# ---------------------------------------------------------------------------
# 7.4 Inicialización de tabla de selección
# ---------------------------------------------------------------------------

table_selected <- data.frame(
  Bioma = character(),
  Month = integer(),
  n_training_rows = integer(),
  n_valid_years_ge30 = integer(),
  Predictors_initial = character(),
  n_predictors_initial = integer(),
  Predictors_after_autocorrelation = character(),
  Formula_after_autocorrelation = character(),
  n_predictors_after_autocorrelation = integer(),
  Predictors_after_RFE = character(),
  Formula_after_RFE = character(),
  n_predictors_after_RFE = integer(),
  stringsAsFactors = FALSE
)

candidate_predictors <- c(
  "f5",
  "count_ActiveFire",
  "prec",
  "temp",
  "FRPsum",
  "FRPmedian",
  "NDVI",
  "FWI",
  "wind",
  "lat",
  "lon",
  "cloud",
  "vpd",
  "soil"
)

response <- "f3"

# ---------------------------------------------------------------------------
# 7.5 Entrenamiento mensual: SOLO selección de predictores
# ---------------------------------------------------------------------------

for (mes in 1:12) {
  
  cat("\nSeleccionando predictores para mes:", mes, "\n")
  
  mes_indices_central <- which(lubridate::month(dates_common) == mes)
  
  mes_indices_window <- sort(unique(c(
    mes_indices_central,
    mes_indices_central - 1,
    mes_indices_central + 1
  )))
  
  mes_indices_window <- mes_indices_window[
    mes_indices_window >= 1 &
      mes_indices_window <= length(dates_common)
  ]
  
  monthly_data <- list()
  
  for (tt in mes_indices_window) {
    
    f3_layer <- arr_common$s3[, , tt]
    f5_layer <- arr_common$f51[, , tt]
    af_layer <- arr_common$af[, , tt]
    
    idx <- which(
      is.finite(f3_layer) &
        is.finite(f5_layer) &
        is.finite(af_layer)
    )
    
    if (length(idx) > 0) {
      
      current_year <- lubridate::year(dates_common[tt])
      key <- as.character(current_year)
      
      df_new <- build_prediction_df(
        arr = arr_common,
        tt = tt,
        idx = idx,
        include_response = TRUE
      )
      
      df_new$year <- current_year
      df_new$month <- lubridate::month(dates_common[tt])
      
      if (!is.null(monthly_data[[key]])) {
        monthly_data[[key]] <- rbind(monthly_data[[key]], df_new)
      } else {
        monthly_data[[key]] <- df_new
      }
    }
  }
  
  if (length(monthly_data) == 0) {
    
    cat("Mes", mes, ": ventana sin datos. No se guarda selección.\n")
    next
  }
  
  df_full <- do.call(rbind, monthly_data)
  
  df_full <- df_full[
    complete.cases(df_full[, c(response, candidate_predictors), drop = FALSE]),
  ]
  
  if (nrow(df_full) == 0) {
    cat("Mes", mes, ": sin casos completos. No se guarda selección.\n")
    next
  }
  
  year_counts <- table(df_full$year)
  n_valid_years_ge30 <- sum(year_counts >= 30)
  
  if (n_valid_years_ge30 < 2) {
    cat("Mes", mes, ": datos insuficientes por año. No se guarda selección.\n")
    next
  }
  
  # -------------------------------------------------------------------------
  # 1) Predictores iniciales válidos: varianza > 0
  # -------------------------------------------------------------------------
  
  valid_predictors <- candidate_predictors
  
  var_check <- sapply(
    df_full[, valid_predictors, drop = FALSE],
    function(x) stats::var(x, na.rm = TRUE)
  )
  
  valid_predictors <- valid_predictors[
    is.finite(var_check) & var_check > 0
  ]
  
  if (length(valid_predictors) == 0) {
    cat("Mes", mes, ": ningún predictor con varianza > 0.\n")
    next
  }
  
  # -------------------------------------------------------------------------
  # 2) Eliminación de autocorrelación por clustering
  # -------------------------------------------------------------------------
  
  if (length(valid_predictors) == 1) {
    
    selected_predictors <- valid_predictors
    
  } else {
    
    cor_matrix <- suppressWarnings(stats::cor(
      df_full[, valid_predictors, drop = FALSE],
      method = "spearman",
      use = "complete.obs"
    ))
    
    cor_matrix[!is.finite(cor_matrix)] <- 0
    diag(cor_matrix) <- 1
    
    dist_matrix <- stats::as.dist(1 - abs(cor_matrix))
    hc <- stats::hclust(dist_matrix, method = "complete")
    groups <- stats::cutree(hc, h = 0.3)
    
    selected_predictors <- c()
    
    for (group in unique(groups)) {
      
      group_vars <- names(groups[groups == group])
      
      corrs <- sapply(group_vars, function(x) {
        cc <- suppressWarnings(stats::cor(
          df_full[[x]],
          df_full[[response]],
          method = "spearman",
          use = "complete.obs"
        ))
        
        if (!is.finite(cc)) cc <- -Inf
        abs(cc)
      })
      
      if (all(!is.finite(corrs))) {
        best_var <- group_vars[1]
      } else {
        best_var <- group_vars[which.max(corrs)]
      }
      
      selected_predictors <- c(selected_predictors, best_var)
    }
  }
  
  selected_predictors <- unique(selected_predictors)
  
  if (length(selected_predictors) == 0) {
    selected_predictors <- valid_predictors
  }
  
  # -------------------------------------------------------------------------
  # 3) RFE sobre los predictores tras autocorrelación
  # -------------------------------------------------------------------------
  
  rfe_ctrl <- caret::rfeControl(
    functions = caret::rfFuncs,
    method = "cv",
    number = 5,
    saveDetails = TRUE,
    verbose = FALSE
  )
  
  set.seed(SEED + biome_id * 1000 + mes)
  
  rfe_fit <- tryCatch(
    {
      caret::rfe(
        x = df_full[, selected_predictors, drop = FALSE],
        y = df_full[[response]],
        sizes = seq(1, length(selected_predictors), by = 1),
        rfeControl = rfe_ctrl,
        ntree = 250,
        mtry = max(1, round(length(selected_predictors) / 3))
      )
    },
    error = function(e) {
      cat("RFE falló en mes", mes, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  
  if (!is.null(rfe_fit)) {
    best_preds <- caret::predictors(rfe_fit)
  } else {
    best_preds <- selected_predictors
  }
  
  best_preds <- unique(best_preds)
  
  if (length(best_preds) == 0) {
    best_preds <- selected_predictors
  }
  
  # -------------------------------------------------------------------------
  # 4) Guardar las dos etapas correctamente
  # -------------------------------------------------------------------------
  
  formula_after_autocorr_txt <- paste(
    response,
    "~",
    paste(selected_predictors, collapse = " + ")
  )
  
  formula_after_RFE_txt <- paste(
    response,
    "~",
    paste(best_preds, collapse = " + ")
  )
  
  table_selected <- rbind(
    table_selected,
    data.frame(
      Bioma = safe_biome_name,
      Month = mes,
      n_training_rows = nrow(df_full),
      n_valid_years_ge30 = n_valid_years_ge30,
      Predictors_initial = paste(valid_predictors, collapse = " + "),
      n_predictors_initial = length(valid_predictors),
      Predictors_after_autocorrelation = paste(selected_predictors, collapse = " + "),
      Formula_after_autocorrelation = formula_after_autocorr_txt,
      n_predictors_after_autocorrelation = length(selected_predictors),
      Predictors_after_RFE = paste(best_preds, collapse = " + "),
      Formula_after_RFE = formula_after_RFE_txt,
      n_predictors_after_RFE = length(best_preds),
      stringsAsFactors = FALSE
    )
  )
  
  cat("Mes", mes, "\n")
  cat("  Iniciales:", length(valid_predictors), "\n")
  cat("  Tras autocorrelación:", length(selected_predictors), "\n")
  cat("  Tras RFE:", length(best_preds), "\n")
  
  rm(df_full, monthly_data, rfe_fit)
  gc()
}

# ---------------------------------------------------------------------------
# 7.6 Guardar CSV del bioma
# ---------------------------------------------------------------------------

out_csv <- file.path(
  output_dir_csv,
  paste0("SelectedPredictors_", safe_biome_name, "_COMMON.csv")
)

write.csv(
  table_selected,
  file = out_csv,
  row.names = FALSE
)

cat("\nCSV guardado en:\n")
cat(out_csv, "\n")

n_months_selected <- nrow(table_selected)

t_end <- Sys.time()

cat("============================================================\n")
cat("Bioma terminado:", bioma, "\n")
cat("Meses con selección guardada:", n_months_selected, "/ 12\n")
cat("Fin:", as.character(t_end), "\n")
cat("Duración segundos:", as.numeric(difftime(t_end, t_start, units = "secs")), "\n")
cat("============================================================\n")

rm(arr, arr_common)
gc()

list(
  Biome = safe_biome_name,
  Status = "OK",
  n_months_selected = n_months_selected,
  seconds = as.numeric(difftime(t_end, t_start, units = "secs")),
  csv = out_csv
)
}