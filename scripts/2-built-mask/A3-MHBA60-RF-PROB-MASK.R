rm(list = ls())
graphics.off()
invisible(gc())

# ===============================
# LIBRERÍAS
# ===============================
suppressPackageStartupMessages({
  library(fields)
  library(maps)
})

# ===============================
# CONFIGURACIÓN GENERAL
# ===============================
Modelo <- "MRBA60-2003-2024-V1"

data_dir <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"
results_root <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask"
output_dir <- file.path(results_root, Modelo, "MASK_FIRE_PREPROC")
maps_dir <- file.path(output_dir, "MAPS_ANNUAL_STATUS")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "RData"), recursive = TRUE, showWarnings = FALSE)
dir.create(maps_dir, recursive = TRUE, showWarnings = FALSE)

# Solo para título del mapa
conf_thresholds <- 30
angle_levels_deg <- 30

# ===============================
# RUTAS DE ENTRADA
# ===============================
file_firecci51 <- file.path(data_dir, "FireCCI51_2003_2024_0.25degree.RData")
file_fireccis311 <- file.path(data_dir, "FireCCIS311_2019_2024_0.25degree.RData")
file_lat <- file.path(data_dir, "latitude.RData")
file_lon <- file.path(data_dir, "longitude.RData")

file_mask_filtered <- file.path(
  output_dir,
  "RData",
  "MASKS_YOUDEN_ONLYFRP_MASK4_and_FRPfiltered.RData"
)

# ===============================
# COMPROBACIÓN DE ARCHIVOS
# ===============================
files_to_check <- c(
  file_firecci51,
  file_fireccis311,
  file_lat,
  file_lon,
  file_mask_filtered
)

missing_files <- files_to_check[!file.exists(files_to_check)]

if (length(missing_files) > 0) {
  stop(
    "Faltan los siguientes archivos:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ===============================
# HELPER PARA CARGAR RDATA
# ===============================
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_rdata.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/to_lon_lat_time2.R")

# ===============================
# FECHAS
# ===============================
dates_full <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_s311 <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
years <- 2003:2024

# ===============================
# CARGA DE DATOS
# ===============================
lon <- read_rdata(file_lon)
lat <- read_rdata(file_lat)

array_fcci51_raw <- read_rdata(file_firecci51) / 1e6
array_s311_raw <- read_rdata(file_fireccis311) / 1e6

count_ActiveFire_tot_filtered_raw <- read_rdata(
  file_mask_filtered,
  "count_ActiveFire_tot_filtered"
)

# ===============================
# NORMALIZAR ORDEN A [lon, lat, time]
# ===============================
array_fcci51 <- to_lon_lat_time2(array_fcci51_raw, lon, lat)
array_s311_in <- to_lon_lat_time2(array_s311_raw, lon, lat)
count_ActiveFire_tot_filtered <- to_lon_lat_time2(
  count_ActiveFire_tot_filtered_raw,
  lon,
  lat
)

rm(array_fcci51_raw, array_s311_raw, count_ActiveFire_tot_filtered_raw)
invisible(gc())


# ===============================
# EXPANDIR FireCCIS311 A 2003-2024
# ===============================
array_s311 <- array(NA_real_, dim = dim(array_fcci51))

ind_s311 <- match(dates_s311, dates_full)

if (any(is.na(ind_s311))) {
  stop("Algunas fechas de FireCCIS311 no coinciden con dates_full.")
}

array_s311[, , ind_s311] <- array_s311_in

rm(array_s311_in)
invisible(gc())

# ===============================
# PREPARAR ARRAYS
# ===============================
array_fcci51[array_fcci51 == 0] <- NA
array_s311[array_s311 == 0] <- NA
count_ActiveFire_tot_filtered[count_ActiveFire_tot_filtered == 0] <- NA

count_ActiveFire_tot <- count_ActiveFire_tot_filtered

# ===============================
# DIMENSIONES
# ===============================
n_months <- dim(array_fcci51)[3]
nrows <- dim(array_fcci51)[1]   # lon
ncols <- dim(array_fcci51)[2]   # lat

# ===============================
# 1) CONSTRUIR FireMask_AF3030F MENSUAL
# ===============================
# Regla corregida:
# FireMask_AF3030F = 1 si CUALQUIERA de los tres datasets tiene valor > 0:
#   - Active Fire filtrado > 0
#   - FireCCI51 > 0
#   - FireCCIS311 > 0
# FireMask_AF3030F = 0 si ningún dataset tiene detección.

FireMask_AF3030F <- array(0L, dim = c(nrows, ncols, n_months))

for (tt in seq_len(n_months)) {
  
  layer_fcci51 <- array_fcci51[, , tt]
  layer_s311 <- array_s311[, , tt]
  layer_mask <- count_ActiveFire_tot[, , tt]
  
  s_present <- !is.na(layer_mask) & layer_mask > 0
  f_present <- !is.na(layer_fcci51) & layer_fcci51 > 0
  f_present3 <- !is.na(layer_s311) & layer_s311 > 0
  
  FireMask_AF3030F[, , tt] <- ifelse(
    s_present | f_present | f_present3,
    1L,
    0L
  )
  
  if (tt %% 12 == 0) {
    cat("Máscara mensual construida hasta:", format(dates_full[tt], "%Y-%m"), "\n")
  }
}

# ===============================
# COMPROBACIÓN DE LA MÁSCARA
# ===============================
any_dataset_present <- 
  (!is.na(array_fcci51) & array_fcci51 > 0) |
  (!is.na(array_s311) & array_s311 > 0) |
  (!is.na(count_ActiveFire_tot) & count_ActiveFire_tot > 0)

check_diff <- FireMask_AF3030F != as.integer(any_dataset_present)

n_diff <- sum(check_diff, na.rm = TRUE)

cat("Número de píxeles/mes distintos entre FireMask_AF3030F y presencia en cualquier dataset:", n_diff, "\n")

if (n_diff != 0) {
  warning("Hay diferencias entre FireMask_AF3030F y la presencia combinada de los datasets.")
} else {
  cat("Comprobación correcta: FireMask_AF3030F = 1 cuando cualquier dataset es > 0.\n")
}

rm(any_dataset_present, check_diff)
invisible(gc())

# ===============================
# 2) GUARDAR FireMask_AF3030F
# ===============================
save(
  FireMask_AF3030F,
  dates_full,
  file = file.path(output_dir, "RData", "FireMask_AF3030F.RData")
)

cat(
  "FireMask_AF3030F guardado en:",
  file.path(output_dir, "RData", "FireMask_AF3030F.RData"),
  "\n"
)

# ===============================
# 3) MAPAS ANUALES DE status_matrix2
# ===============================
for (conf in conf_thresholds) {
  for (angle_deg in angle_levels_deg) {
    for (year in years) {
      
      cat("Procesando mapa anual:", year, "\n")
      
      year_idx <- which(format(dates_full, "%Y") == as.character(year))
      
      if (length(year_idx) != 12) {
        warning("El año ", year, " no tiene 12 meses completos. Se omite.")
        next
      }
      
      total_fcci51 <- apply(
        array_fcci51[, , year_idx, drop = FALSE],
        c(1, 2),
        sum,
        na.rm = TRUE
      )
      
      total_s311 <- apply(
        array_s311[, , year_idx, drop = FALSE],
        c(1, 2),
        sum,
        na.rm = TRUE
      )
      
      total_mask <- apply(
        count_ActiveFire_tot[, , year_idx, drop = FALSE],
        c(1, 2),
        sum,
        na.rm = TRUE
      )
      
      total_fcci51[total_fcci51 == 0] <- NA
      total_s311[total_s311 == 0] <- NA
      total_mask[total_mask == 0] <- NA
      
      # Presencia anual
      s_present <- !is.na(total_mask) & total_mask > 0
      f_present <- !is.na(total_fcci51) & total_fcci51 > 0
      f_present3 <- !is.na(total_s311) & total_s311 > 0
      
      # Clasificación anual
      status_matrix2 <- matrix(NA_real_, nrow = nrows, ncol = ncols)
      
      status_matrix2[s_present  & !f_present & !f_present3] <- 0  # Solo Active Fire
      status_matrix2[!s_present & !f_present &  f_present3] <- 1  # Solo FireCCIS311
      status_matrix2[!s_present &  f_present & !f_present3] <- 2  # Solo FireCCI51
      status_matrix2[s_present  & !f_present &  f_present3] <- 3  # Active Fire + FireCCIS311
      status_matrix2[!s_present &  f_present &  f_present3] <- 4  # FireCCI51 + FireCCIS311
      status_matrix2[s_present  &  f_present & !f_present3] <- 5  # Active Fire + FireCCI51
      status_matrix2[s_present  &  f_present &  f_present3] <- 6  # Los tres
      
      # Conteos
      count_only_active <- sum(status_matrix2 == 0, na.rm = TRUE)
      count_only_s3 <- sum(status_matrix2 == 1, na.rm = TRUE)
      count_only_fcci51 <- sum(status_matrix2 == 2, na.rm = TRUE)
      count_active_s3 <- sum(status_matrix2 == 3, na.rm = TRUE)
      count_fcci51_s3 <- sum(status_matrix2 == 4, na.rm = TRUE)
      count_active_fcci51 <- sum(status_matrix2 == 5, na.rm = TRUE)
      count_all_three <- sum(status_matrix2 == 6, na.rm = TRUE)
      
      total_pixels <- sum(status_matrix2 %in% 0:6, na.rm = TRUE)
      
      if (total_pixels == 0) {
        warning("No hay píxeles con detección en ", year, ". Se omite el mapa.")
        next
      }
      
      perc_only_active <- (count_only_active / total_pixels) * 100
      perc_only_s3 <- (count_only_s3 / total_pixels) * 100
      perc_only_fcci51 <- (count_only_fcci51 / total_pixels) * 100
      perc_active_s3 <- (count_active_s3 / total_pixels) * 100
      perc_fcci51_s3 <- (count_fcci51_s3 / total_pixels) * 100
      perc_active_fcci51 <- (count_active_fcci51 / total_pixels) * 100
      perc_all_three <- (count_all_three / total_pixels) * 100
      
      # ==============================
      # RECORTE DE LATITUDES
      # ==============================
      lat_indices <- which(lat >= -58 & lat <= 85)
      
      if (length(lat_indices) == 0) {
        stop("No hay latitudes dentro del rango [-58, 85].")
      }
      
      lat_filtered <- lat[lat_indices]
      status_matrix_filtered <- status_matrix2[, lat_indices, drop = FALSE]
      
      # Ordenar para que image() reciba ejes crecientes
      lon_order <- order(lon)
      lat_order <- order(lat_filtered)
      
      lon_plot <- lon[lon_order]
      lat_plot <- lat_filtered[lat_order]
      status_matrix_plot <- status_matrix_filtered[lon_order, lat_order, drop = FALSE]
      
      # ==============================
      # COLORES Y CORTES
      # ==============================
      my_breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)
      
      my_colors <- c(
        "blue",    # 0 Solo Active Fire
        "red",     # 1 Solo FireCCIS311
        "green",   # 2 Solo FireCCI51
        "orange",  # 3 Active Fire + FireCCIS311
        "purple",  # 4 FireCCI51 + FireCCIS311
        "cyan",    # 5 Active Fire + FireCCI51
        "yellow"   # 6 Los tres
      )
      
      output_file <- file.path(
        maps_dir,
        paste0("mapa_all_", year, "_WGS84_025_FILTER.jpg")
      )
      
      jpeg(file = output_file, width = 12, height = 8, units = "in", res = 900)
      
      par(mar = c(6.2, 3.5, 3.2, 1.2), xpd = NA)
      
      image(
        x = lon_plot,
        y = lat_plot,
        z = status_matrix_plot,
        col = my_colors,
        breaks = my_breaks,
        axes = FALSE,
        xlab = "",
        ylab = "",
        useRaster = TRUE
      )
      
      axis(1)
      axis(2, las = 1)
      box()
      
      title(
        main = paste0(
          "Fire detection (", year,
          "): Active fires, FireCCI51, FireCCIS311; CL \u2265 ",
          conf, ", \u03B8 \u2264 ", angle_deg, "\u00B0"
        ),
        cex.main = 1.4,
        line = 0.8
      )
      
      # Namespace explícito para evitar conflicto con purrr::map()
      maps::map(
        "world",
        add = TRUE,
        col = "black",
        lwd = 0.5
      )
      
      abline(h = c(-40, 40), col = "gray28", lwd = 1)
      abline(h = 0, col = "gray28", lwd = 1.5)
      
      legend_text <- c(
        paste0("Only Active fires: ", round(perc_only_active, 1), "%"),
        paste0("Only FireCCIS311: ", round(perc_only_s3, 1), "%"),
        paste0("Only FireCCI51: ", round(perc_only_fcci51, 1), "%"),
        paste0("Active fires + FireCCIS311: ", round(perc_active_s3, 1), "%"),
        paste0("FireCCI51 + FireCCIS311: ", round(perc_fcci51_s3, 1), "%"),
        paste0("Active fires + FireCCI51: ", round(perc_active_fcci51, 1), "%"),
        paste0("Three products: ", round(perc_all_three, 1), "%")
      )
      
      legend(
        x = "bottom",
        legend = legend_text,
        fill = my_colors,
        bty = "n",
        cex = 0.95,
        ncol = 3,
        inset = c(0, -0.18)
      )
      
      par(xpd = FALSE)
      dev.off()
      
      cat("Mapa guardado en:", output_file, "\n")
      
      rm(
        total_fcci51,
        total_s311,
        total_mask,
        s_present,
        f_present,
        f_present3,
        status_matrix2,
        status_matrix_filtered,
        status_matrix_plot
      )
      
      invisible(gc())
    }
  }
}
