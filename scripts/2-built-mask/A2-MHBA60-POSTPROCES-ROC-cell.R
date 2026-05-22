# ===============================
# CONFIGURACIÓN PREVIA
# ===============================
rm(list = ls())
graphics.off()
invisible(gc())

library(pROC)
library(dplyr)
library(fields)
library(sf)
library(lubridate)
library(ggplot2)
library(ggpmisc)
library(viridis)
library(ggpointdensity)
library(maps)

# ===============================
# CONFIGURACIÓN GENERAL
# ===============================
Modelo <- "MRBA60-2003-2024-V1"

# Ruta base nueva
data_dir <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

# Ruta shapefile biomas
file_biomes <- "/mnt/disco6tb/MRBA60/data/A1_RAW/MBC/continental-biomes_dinerstein_V10.shp"

# Active Fire disponible en A3_ADJ:
# "FRPsum", "FRPmean" o "AFcount"
active_fire_metric <- "AFcount"

# Salidas
results_root <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask"
output_dir <- file.path(results_root, Modelo, "MASK_FIRE_PREPROC")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "RData"), recursive = TRUE, showWarnings = FALSE)

# Activa esto solo si quieres ver mapas de depuración
plot_debug_masks <- FALSE

# ===============================
# RUTAS DE ENTRADA
# ===============================
file_firecci51 <- file.path(data_dir, "FireCCI51_2003_2024_0.25degree.RData")
file_fireccis311 <- file.path(data_dir, "FireCCIS311_2019_2024_0.25degree.RData")
file_lon <- file.path(data_dir, "longitude.RData")
file_lat <- file.path(data_dir, "latitude.RData")
file_active_fire <- file.path(
  data_dir,
  paste0("MODIS-", active_fire_metric, "_conf30_angle30-200301-202412-025.RData")
)

file_prob_full <- file.path(
  results_root,
  Modelo,
  "RData",
  paste0("BA_", Modelo, "global_BA_PROB_FireHarmonized_Full.RData")
)

# ===============================
# COMPROBACIÓN DE ARCHIVOS
# ===============================
files_to_check <- c(
  file_firecci51,
  file_fireccis311,
  file_lon,
  file_lat,
  file_active_fire,
  file_biomes,
  file_prob_full
)

missing_files <- files_to_check[!file.exists(files_to_check)]
if (length(missing_files) > 0) {
  stop(
    "Faltan los siguientes archivos de entrada:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ===============================
# HELPER PARA CARGAR .RData
# ===============================
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_rdata.R")
sanitize_biome_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("[^[:alnum:]_]", "", x)
  x
}

source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/to_lon_lat_time.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/apply_mask_3d.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/build_biome_mask_info.R")

# ============================================================================
# FECHAS: PERÍODO COMPLETO Y PERÍODO COMÚN
# ============================================================================
dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

# ============================================================================
# LEER BIOMAS
# ============================================================================
biomas_shp <- st_read(file_biomes, quiet = TRUE)
biomas_shp <- st_transform(biomas_shp, crs = 4326)

# ============================================================================
# LEER LONGITUD Y LATITUD
# ============================================================================
lon <- read_rdata(file_lon)
lat <- read_rdata(file_lat)

# ============================================================================
# LEER DATOS RAW
# ============================================================================
# FireCCI51 (2003-2024)
BA_Fire51_tot_raw <- read_rdata(file_firecci51) / 1e6
BA_Fire51_tot_raw[BA_Fire51_tot_raw == 0] <- NA
invisible(gc())

# FireCCIS311 (2019-2024)
BA_FireS3_raw <- read_rdata(file_fireccis311) / 1e6
BA_FireS3_raw[BA_FireS3_raw == 0] <- NA
invisible(gc())

# Active Fire
count_ActiveFire_tot_raw <- read_rdata(file_active_fire)
count_ActiveFire_tot_raw[count_ActiveFire_tot_raw == 0] <- NA
invisible(gc())

# Probabilidad armonizada
e_prob <- new.env(parent = emptyenv())
nm_prob <- load(file_prob_full, envir = e_prob)

if ("global_BA_FireHarmonized_prob_full" %in% nm_prob) {
  global_BA_FireHarmonized_prob_full_raw <- e_prob$global_BA_FireHarmonized_prob_full
} else if (length(nm_prob) == 1) {
  global_BA_FireHarmonized_prob_full_raw <- e_prob[[nm_prob[1]]]
} else {
  stop(
    "El archivo de probabilidad contiene varios objetos y no aparece ",
    "'global_BA_FireHarmonized_prob_full'. Objetos: ",
    paste(nm_prob, collapse = ", ")
  )
}
rm(e_prob, nm_prob)
invisible(gc())


# ============================================================================
# DETECTAR ORIENTACIÓN ESPACIAL
# ============================================================================
dims_f5 <- dim(BA_Fire51_tot_raw)

cat("Dimensiones FireCCI51:", paste(dims_f5, collapse = " x "), "\n")
cat("Longitud lon:", length(lon), "\n")
cat("Longitud lat:", length(lat), "\n")

if (dims_f5[1] == length(lat) && dims_f5[2] == length(lon)) {
  grid_order <- "lat_lon_time"
} else if (dims_f5[1] == length(lon) && dims_f5[2] == length(lat)) {
  grid_order <- "lon_lat_time"
} else {
  stop(
    "Las dimensiones espaciales no coinciden con lon/lat.\n",
    "dim(BA_Fire51_tot_raw) = ", paste(dims_f5, collapse = " x "), "\n",
    "length(lon) = ", length(lon), "\n",
    "length(lat) = ", length(lat)
  )
}

cat("Orden espacial detectado:", grid_order, "\n")

# ============================================================================
# ESTANDARIZAR TODO A [lon, lat, time]
# ============================================================================
BA_Fire51_tot <- to_lon_lat_time(BA_Fire51_tot_raw, grid_order)
BA_FireS3     <- to_lon_lat_time(BA_FireS3_raw, grid_order)
count_ActiveFire_tot <- to_lon_lat_time(count_ActiveFire_tot_raw, grid_order)
global_BA_FireHarmonized_prob_full <- to_lon_lat_time(global_BA_FireHarmonized_prob_full_raw, grid_order)

rm(BA_Fire51_tot_raw, BA_FireS3_raw, count_ActiveFire_tot_raw, global_BA_FireHarmonized_prob_full_raw)
invisible(gc())

# ============================================================================
# EXPANDIR FireCCIS311 a 2003-2024
# ============================================================================
nlon <- dim(BA_Fire51_tot)[1]
nlat <- dim(BA_Fire51_tot)[2]
ntime_full <- dim(BA_Fire51_tot)[3]

BA_FireS3_tot <- array(NA_real_, dim = c(nlon, nlat, ntime_full))
BA_FireS3_tot[, , ind_common] <- BA_FireS3
rm(BA_FireS3)
invisible(gc())

# ============================================================================
# ARRAY BASE A FILTRAR
# ============================================================================
global_BA_FireHarmonized_full <- count_ActiveFire_tot

# ============================================================================
# MÁSCARA PARA CELDAS VÁLIDAS
# ============================================================================
s_present  <- !is.na(count_ActiveFire_tot) & count_ActiveFire_tot > 0
f5_present <- !is.na(BA_Fire51_tot)        & BA_Fire51_tot > 0
f3_present <- !is.na(BA_FireS3_tot)        & BA_FireS3_tot > 0

status_matrix2_tot <- array(NA_real_, dim = dim(BA_Fire51_tot))
status_matrix2_tot[s_present | f5_present | f3_present] <- 1

rm(s_present, f5_present, f3_present)
invisible(gc())

# ============================================================================
# COMPLETAR CEROS EN CELDAS VÁLIDAS
# ============================================================================
global_BA_FireHarmonized_prob_full[
  status_matrix2_tot == 1 & is.na(global_BA_FireHarmonized_prob_full)
] <- 0

BA_Fire51_tot[
  status_matrix2_tot == 1 & is.na(BA_Fire51_tot)
] <- 0

BA_FireS3_tot[
  status_matrix2_tot == 1 & is.na(BA_FireS3_tot)
] <- 0

count_ActiveFire_tot[
  status_matrix2_tot == 1 & is.na(count_ActiveFire_tot)
] <- 0

global_BA_FireHarmonized_full[
  status_matrix2_tot == 1 & is.na(global_BA_FireHarmonized_full)
] <- 0

# ============================================================================
# MALLADO [lon, lat]
# ============================================================================
lon_mat <- matrix(
  rep(lon, times = length(lat)),
  nrow = length(lon),
  ncol = length(lat),
  byrow = FALSE
)

lat_mat <- matrix(
  rep(lat, each = length(lon)),
  nrow = length(lon),
  ncol = length(lat),
  byrow = FALSE
)

lon_range <- lon_mat[, 1]
lat_range <- lat_mat[1, ]

# ========================
# INICIALIZACIÓN
# ========================
global_BA_FireHarmonized_full_filtered <- global_BA_FireHarmonized_full

roc_results_table <- data.frame(
  Bioma = character(),
  Bioma_safe = character(),
  Mes = integer(),
  Threshold = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  Youden = numeric(),
  AUC = numeric(),
  Pvalue = numeric(),
  Pvalue_FDR = numeric(),
  Threshold_use = numeric(),
  N_valid = integer(),
  N_pos = integer(),
  N_neg = integer(),
  stringsAsFactors = FALSE
)

biomas_unique <- sort(unique(na.omit(biomas_shp$cont_bm)))
biome_info_list <- vector("list", length(biomas_unique))
names(biome_info_list) <- biomas_unique

# ============================================================================
# FASE 1: CALCULAR ROC / THRESHOLD POR BIOMA Y MES
# ============================================================================
for (bioma in biomas_unique) {
  cat("\nProcesando bioma (fase ROC):", bioma, "\n")
  
  biome_info <- build_biome_mask_info(
    bioma = bioma,
    biomas_shp = biomas_shp,
    lon_range = lon_range,
    lat_range = lat_range,
    lon_mat = lon_mat,
    lat_mat = lat_mat,
    status_matrix2_tot = status_matrix2_tot
  )
  
  if (is.null(biome_info)) {
    cat("  Sin celdas válidas para este bioma\n")
    next
  }
  
  biome_info_list[[bioma]] <- biome_info
  
  lon_idx <- biome_info$lon_idx
  lat_idx <- biome_info$lat_idx
  mask_matrix <- biome_info$mask_matrix
  
  if (!any(mask_matrix, na.rm = TRUE)) {
    cat("  Máscara vacía para este bioma\n")
    next
  }
  
  harmonized_crop <- global_BA_FireHarmonized_full[lon_idx, lat_idx, ind_common, drop = FALSE]
  prob_crop       <- global_BA_FireHarmonized_prob_full[lon_idx, lat_idx, ind_common, drop = FALSE]
  fire51_crop     <- BA_Fire51_tot[lon_idx, lat_idx, ind_common, drop = FALSE]
  fireS3_crop     <- BA_FireS3_tot[lon_idx, lat_idx, ind_common, drop = FALSE]
  
  harmonized_crop <- apply_mask_3d(harmonized_crop, mask_matrix)
  prob_crop       <- apply_mask_3d(prob_crop, mask_matrix)
  fire51_crop     <- apply_mask_3d(fire51_crop, mask_matrix)
  fireS3_crop     <- apply_mask_3d(fireS3_crop, mask_matrix)
  
  if (plot_debug_masks) {
    image.plot(
      lon_range[lon_idx],
      lat_range[lat_idx],
      mask_matrix,
      main = paste("Máscara:", bioma),
      xlab = "Lon",
      ylab = "Lat"
    )
    plot(st_geometry(biome_info$bioma_sel), add = TRUE, border = "black", lwd = 1.5)
  }
  
  for (mes in 1:12) {
    cat("  Procesando mes:", mes, "\n")
    
    mes_indices <- which(month(dates_common) == mes)
    if (length(mes_indices) == 0) next
    
    harmonized_mes <- harmonized_crop[, , mes_indices, drop = FALSE]
    prob_mes       <- prob_crop[, , mes_indices, drop = FALSE]
    fire51_mes     <- fire51_crop[, , mes_indices, drop = FALSE]
    fireS3_mes     <- fireS3_crop[, , mes_indices, drop = FALSE]
    
    prob_vec       <- as.vector(prob_mes)
    harmonized_vec <- as.vector(harmonized_mes)
    fire51_vec     <- as.vector(fire51_mes)
    fireS3_vec     <- as.vector(fireS3_mes)
    
    valid_idx <- which(!is.na(prob_vec) & !is.na(harmonized_vec))
    
    if (length(valid_idx) < 30) {
      cat("    Insuficientes datos válidos en mes", mes, "\n")
      next
    }
    
    has_fire51 <- !is.na(fire51_vec[valid_idx]) & (fire51_vec[valid_idx] > 0)
    has_fireS3 <- !is.na(fireS3_vec[valid_idx]) & (fireS3_vec[valid_idx] > 0)
    
    labels <- ifelse(has_fire51 | has_fireS3, 1L, 0L)
    
    n_pos <- sum(labels == 1L, na.rm = TRUE)
    n_neg <- sum(labels == 0L, na.rm = TRUE)
    
    if (n_pos == 0 || n_neg == 0) {
      cat("    No hay suficiente variabilidad en las etiquetas\n")
      next
    }
    
    roc_obj <- tryCatch(
      roc(labels, prob_vec[valid_idx], quiet = TRUE),
      error = function(e) NULL
    )
    
    if (is.null(roc_obj)) {
      cat("    Error calculando ROC\n")
      next
    }
    
    auc_value <- as.numeric(auc(roc_obj))
    
    ci_auc <- tryCatch(
      ci.auc(roc_obj),
      error = function(e) c(NA_real_, NA_real_, NA_real_)
    )
    
    if (all(is.finite(ci_auc))) {
      se_auc <- (ci_auc[3] - ci_auc[1]) / (2 * 1.96)
      if (is.finite(se_auc) && se_auc > 0) {
        z_value <- (auc_value - 0.5) / se_auc
        pvalue <- 2 * pnorm(-abs(z_value))
      } else {
        pvalue <- NA_real_
      }
    } else {
      pvalue <- NA_real_
    }
    
    coords_best <- coords(
      roc_obj,
      x = "best",
      best.method = "youden",
      ret = c("threshold", "sensitivity", "specificity"),
      transpose = FALSE
    )
    
    if (is.list(coords_best)) coords_best <- unlist(coords_best)
    
    threshold_opt <- as.numeric(coords_best["threshold"])
    sens_opt      <- as.numeric(coords_best["sensitivity"])
    spec_opt      <- as.numeric(coords_best["specificity"])
    youden_opt    <- sens_opt + spec_opt - 1
    
    if (length(threshold_opt) == 0 || is.na(threshold_opt) || !is.finite(threshold_opt)) {
      cat("    No se encontró threshold óptimo\n")
      next
    }
    
    new_result <- data.frame(
      Bioma = bioma,
      Bioma_safe = biome_info$safe_biome_name,
      Mes = mes,
      Threshold = threshold_opt,
      Sensitivity = sens_opt,
      Specificity = spec_opt,
      Youden = youden_opt,
      AUC = auc_value,
      Pvalue = pvalue,
      Pvalue_FDR = NA_real_,
      Threshold_use = NA_real_,
      N_valid = length(valid_idx),
      N_pos = n_pos,
      N_neg = n_neg,
      stringsAsFactors = FALSE
    )
    
    roc_results_table <- rbind(roc_results_table, new_result)
  }
}

# ============================================================================
# APLICAR FDR GLOBAL
# ============================================================================
if (nrow(roc_results_table) == 0) {
  warning("No se generaron resultados ROC. El producto filtrado será idéntico al original.")
} else {
  roc_results_table$Pvalue_FDR <- p.adjust(roc_results_table$Pvalue, method = "fdr")
  roc_results_table$Threshold_use <- ifelse(
    !is.na(roc_results_table$Pvalue_FDR) & roc_results_table$Pvalue_FDR < 0.05,
    roc_results_table$Threshold,
    0
  )
}

# ============================================================================
# FASE 2: APLICAR FILTRO POR BIOMA Y MES A TODO EL PERIODO
# ============================================================================
if (nrow(roc_results_table) > 0) {
  for (bioma in biomas_unique) {
    biome_info <- biome_info_list[[bioma]]
    
    if (is.null(biome_info)) next
    
    cat("\nAplicando filtro al bioma:", bioma, "\n")
    
    lon_idx <- biome_info$lon_idx
    lat_idx <- biome_info$lat_idx
    mask_matrix <- biome_info$mask_matrix
    
    biome_thresholds <- roc_results_table %>%
      filter(Bioma == bioma)
    
    if (nrow(biome_thresholds) == 0) {
      cat("  Sin thresholds disponibles para este bioma\n")
      next
    }
    
    for (ii in seq_len(nrow(biome_thresholds))) {
      mes <- biome_thresholds$Mes[ii]
      threshold_use <- biome_thresholds$Threshold_use[ii]
      
      cat("  Mes:", mes, " | threshold_use =", threshold_use, "\n")
      
      if (is.na(threshold_use) || !is.finite(threshold_use) || threshold_use <= 0) next
      
      mes_indices_full <- which(month(dates_full) == mes)
      
      for (t in mes_indices_full) {
        harmonized_layer <- global_BA_FireHarmonized_full_filtered[lon_idx, lat_idx, t, drop = TRUE]
        prob_layer       <- global_BA_FireHarmonized_prob_full[lon_idx, lat_idx, t, drop = TRUE]
        fire51_layer     <- BA_Fire51_tot[lon_idx, lat_idx, t, drop = TRUE]
        fireS3_layer     <- BA_FireS3_tot[lon_idx, lat_idx, t, drop = TRUE]
        
        idx_valid <- which(!is.na(prob_layer) & mask_matrix)
        if (length(idx_valid) == 0) next
        
        vprob <- prob_layer[idx_valid]
        
        has_fire51 <- !is.na(fire51_layer[idx_valid]) & (fire51_layer[idx_valid] > 0)
        has_fireS3 <- !is.na(fireS3_layer[idx_valid]) & (fireS3_layer[idx_valid] > 0)
        
        year_current <- year(dates_full[t])
        
        if (year_current <= 2018) {
          has_other_detection <- has_fire51
        } else {
          has_other_detection <- has_fire51 | has_fireS3
        }
        
        has_other_detection[is.na(has_other_detection)] <- FALSE
        
        low_prob <- vprob < threshold_use
        to_zero_logical <- low_prob & (!has_other_detection)
        
        idx_false <- idx_valid[to_zero_logical]
        
        if (length(idx_false) > 0) {
          harmonized_layer[idx_false] <- 0
          global_BA_FireHarmonized_full_filtered[lon_idx, lat_idx, t] <- harmonized_layer
        }
      }
    }
  }
}

# ============================================================================
# GUARDAR RESULTADOS
# ============================================================================
count_ActiveFire_tot_filtered <- global_BA_FireHarmonized_full_filtered

save(
  count_ActiveFire_tot_filtered,
  file = file.path(output_dir, "RData", "MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData")
)

write.csv(
  roc_results_table,
  file.path(output_dir, "csv", "ROC_Thresholds_ByBioma_Month_MASK.csv"),
  row.names = FALSE
)

save(
  roc_results_table,
  file = file.path(output_dir, "RData", "ROC_Thresholds_ByBioma_Month_MASK.RData")
)
