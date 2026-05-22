rm(list = ls())
graphics.off()
invisible(gc())

# ===============================
# LIBRERÍAS
# ===============================
library(fields)
library(maps)

# ===============================
# CONFIGURACIÓN GENERAL
# ===============================
Modelo <- "MRBA60-2003-2024-V1"

data_dir <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"
results_root <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask"

output_dir <- file.path(
  results_root,
  Modelo,
  "MASK_FIRE_PREPROC",
  "plot_mask_NoFilter"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

active_fire_metric <- "AFcount"
conf_thresholds <- 30
angle_levels_deg <- 30

# ===============================
# RUTAS DE ENTRADA
# ===============================
file_firecci51 <- file.path(data_dir, "FireCCI51_2003_2024_0.25degree.RData")
file_fireccis311 <- file.path(data_dir, "FireCCIS311_2019_2024_0.25degree.RData")
file_lat <- file.path(data_dir, "latitude.RData")
file_lon <- file.path(data_dir, "longitude.RData")
file_active_fire <- file.path(
  data_dir,
  paste0("MODIS-", active_fire_metric, "_conf30_angle30-200301-202412-025.RData")
)

# ===============================
# COMPROBACIÓN DE ARCHIVOS
# ===============================
files_to_check <- c(
  file_firecci51,
  file_fireccis311,
  file_lat,
  file_lon,
  file_active_fire
)

source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/missing_files.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/read_rdata.R")
source("/mnt/disco6tb/Dropbox/UAH/FireCCI60/paper/review_v1/scripts/0_functions/to_lon_lat_time2.R")

# ===============================
# FECHAS
# ===============================
dates_full <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_s311 <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")

years <- 2003:2024
# Para mantener exactamente el periodo antiguo:
# years <- 2003:2022

# ===============================
# CARGA DE DATOS
# ===============================
lon <- read_rdata(file_lon)
lat <- read_rdata(file_lat)

array_fcci51_raw <- read_rdata(file_firecci51) / 1e6
array_s311_raw   <- read_rdata(file_fireccis311) / 1e6
count_ActiveFire_tot_raw <- read_rdata(file_active_fire)

count_ActiveFire_tot_raw[count_ActiveFire_tot_raw == 0] <- NA

# ===============================
# NORMALIZAR ORDEN A [lon, lat, time]
# ===============================
array_fcci51 <- to_lon_lat_time2(array_fcci51_raw, lon, lat)
array_s311_in <- to_lon_lat_time2(array_s311_raw, lon, lat)
count_ActiveFire_tot <- to_lon_lat_time2(count_ActiveFire_tot_raw, lon, lat)

rm(array_fcci51_raw, array_s311_raw, count_ActiveFire_tot_raw)
invisible(gc())

# ===============================
# EXPANDIR FireCCIS311 A 2003-2024
# ===============================
array_s311 <- array(NA_real_, dim = dim(array_fcci51))
ind_s311 <- which(dates_full %in% dates_s311)
array_s311[, , ind_s311] <- array_s311_in
rm(array_s311_in)
invisible(gc())

# ===============================
# DIMENSIONES
# ===============================
nrows <- dim(array_fcci51)[1]   # lon
ncols <- dim(array_fcci51)[2]   # lat

# ===============================
# BUCLE PRINCIPAL
# ===============================
for (conf in conf_thresholds) {
  for (angle_deg in angle_levels_deg) {
    for (year in years) {
      
      year_idx <- which(format(dates_full, "%Y") == as.character(year))
      
      if (length(year_idx) != 12) {
        warning("El año ", year, " no tiene 12 meses completos. Se omite.")
        next
      }
      
      total_fcci51 <- apply(array_fcci51[, , year_idx, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      total_s311   <- apply(array_s311[, , year_idx, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      total_mask   <- apply(count_ActiveFire_tot[, , year_idx, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
      
      total_fcci51[total_fcci51 == 0] <- NA
      total_s311[total_s311 == 0] <- NA
      total_mask[total_mask == 0] <- NA
      
      # Presencia anual
      s_present  <- !is.na(total_mask)   & total_mask > 0
      f_present  <- !is.na(total_fcci51) & total_fcci51 > 0
      f_present3 <- !is.na(total_s311)   & total_s311 > 0
      
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
      count_all_three     <- sum(status_matrix2 == 6, na.rm = TRUE)
      count_fcci51_puntos <- sum(status_matrix2 == 5, na.rm = TRUE)
      count_fcci51_s3     <- sum(status_matrix2 == 4, na.rm = TRUE)
      count_active_s3     <- sum(status_matrix2 == 3, na.rm = TRUE)
      count_only_fcci51   <- sum(status_matrix2 == 2, na.rm = TRUE)
      count_only_s3       <- sum(status_matrix2 == 1, na.rm = TRUE)
      count_only_active   <- sum(status_matrix2 == 0, na.rm = TRUE)
      
      total_pixels <- sum(status_matrix2 %in% 0:6, na.rm = TRUE)
      
      if (total_pixels == 0) {
        warning("No hay píxeles con detección en ", year, ". Se omite el mapa.")
        next
      }
      
      # Porcentajes
      perc_all_three     <- (count_all_three / total_pixels) * 100
      perc_fcci51_puntos <- (count_fcci51_puntos / total_pixels) * 100
      perc_fcci51_s3     <- (count_fcci51_s3 / total_pixels) * 100
      perc_active_s3     <- (count_active_s3 / total_pixels) * 100
      perc_only_fcci51   <- (count_only_fcci51 / total_pixels) * 100
      perc_only_s3       <- (count_only_s3 / total_pixels) * 100
      perc_only_active   <- (count_only_active / total_pixels) * 100
      
      # ==============================
      # RECORTE DE LATITUDES
      # ==============================
      lat_indices <- which(lat >= -58 & lat <= 85)
      
      if (length(lat_indices) == 0) {
        stop("No hay latitudes dentro del rango [-58, 85].")
      }
      
      lat_filtered <- lat[lat_indices]
      status_matrix_filtered <- status_matrix2[, lat_indices, drop = FALSE]
      
      # Ordenar ejes para image()
      lon_order <- order(lon)
      lat_order <- order(lat_filtered)
      
      lon_plot <- lon[lon_order]
      lat_plot <- lat_filtered[lat_order]
      status_matrix_plot <- status_matrix_filtered[lon_order, lat_order, drop = FALSE]
      
      # ==============================
      # COLORES Y CORTES
      # ==============================
      my_breaks <- c(-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5)
      my_colors <- c("blue", "red", "green", "orange", "purple", "cyan", "yellow")
      
      output_file <- file.path(
        output_dir,
        paste0("mapa_all_", year, "_WGS84_025_NOFILTER.jpg")
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
        paste0("Active fires + FireCCI51: ", round(perc_fcci51_puntos, 1), "%"),
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
        total_fcci51, total_s311, total_mask,
        s_present, f_present, f_present3,
        status_matrix2, status_matrix_filtered, status_matrix_plot
      )
      invisible(gc())
    }
  }
}

if (dev.cur() > 1) dev.off()

