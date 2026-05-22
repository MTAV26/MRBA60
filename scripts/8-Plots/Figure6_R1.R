# ==========================================================
# SCRIPT COMPLETO CORREGIDO
# - MRBA60, FireCCI51, MCD64A1, GFED5
# - Mapas de media anual de área quemada
# - Mapas de diferencias MRBA60 - producto
# - 0 y NA se dejan SIN COLOR
# - En mapas de media:
#     0 < x < 5       -> primer color
#     5 <= x < 10     -> segundo color
#     ...
#     x >= 750        -> último color
# - En mapas de diferencias:
#     diferencia == 0 -> sin color
#     clases discretas coherentes con breaks
# ==========================================================


# ==========================================================
# Limpieza inicial
# ==========================================================
rm(list = ls())
graphics.off()
gc()


# ==========================================================
# Librerías
# ==========================================================
suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(ncdf4)
  library(sf)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(fields)
  library(viridisLite)
})


# ==========================================================
# Configuración de rutas y modelo
# ==========================================================
Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_base      <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv   <- file.path(output_base, "csv")
output_dir_plot  <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure6_R1/"
output_dir_RData <- file.path(output_base, "RData")

dir.create(output_base,      showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_csv,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_plot,  showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_RData, showWarnings = FALSE, recursive = TRUE)


# ==========================================================
# 1) Parámetros temporales
# ==========================================================
dates_0324   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
years_0324   <- 2003:2024
nmonths_0324 <- length(dates_0324)   # 264 meses

stopifnot(nmonths_0324 == 264)


# ==========================================================
# 2) Funciones auxiliares
# ==========================================================

# ----------------------------------------------------------
# Recortar array 3D [lon, lat, time] a los primeros nkeep meses
# ----------------------------------------------------------
crop_months <- function(arr_3d, nkeep) {
  
  if (length(dim(arr_3d)) != 3) {
    stop("Se esperaba un array 3D [lon, lat, time].")
  }
  
  if (dim(arr_3d)[3] < nkeep) {
    stop(sprintf(
      "El array tiene %d meses y se requieren %d.",
      dim(arr_3d)[3], nkeep
    ))
  }
  
  arr_3d[, , 1:nkeep, drop = FALSE]
}


# ----------------------------------------------------------
# Suma mensual -> suma anual -> media anual
# Devuelve matriz [lon, lat] en km2/año medio
# ----------------------------------------------------------
monthly_to_annual_mean <- function(arr_3d, nkeep) {
  
  arr <- crop_months(arr_3d, nkeep)
  
  nx <- dim(arr)[1]
  ny <- dim(arr)[2]
  
  if (nkeep %% 12 != 0) {
    stop("nkeep debe ser múltiplo de 12.")
  }
  
  nyears <- nkeep / 12
  
  arr4 <- array(arr, dim = c(nx, ny, 12, nyears))
  
  annual_sum <- apply(arr4, c(1, 2, 4), sum, na.rm = TRUE)
  annual_mean <- apply(annual_sum, c(1, 2), mean, na.rm = TRUE)
  
  annual_mean
}


# ----------------------------------------------------------
# Matriz [lon, lat] -> RasterLayer WGS84
# Excluye Antártida: lat > -60
# ----------------------------------------------------------
mat_to_raster_ll <- function(mat, lon, lat, lat_min = -60) {
  
  if (lat[1] > tail(lat, 1)) {
    lat <- rev(lat)
    mat <- mat[, ncol(mat):1]
  }
  
  keep <- which(lat > lat_min)
  
  mat <- mat[, keep, drop = FALSE]
  lat <- lat[keep]
  
  r_mat <- t(mat)
  
  r <- raster::raster(
    r_mat,
    xmn = min(lon),
    xmx = max(lon),
    ymn = min(lat),
    ymx = max(lat),
    crs = sp::CRS("+proj=longlat +datum=WGS84")
  )
  
  raster::flip(r, "y")
}


# ----------------------------------------------------------
# Proyección Robinson
# ----------------------------------------------------------
crs_robinson <- "+proj=robin +lon_0=0 +datum=WGS84 +units=m +no_defs"

to_robin <- function(r_ll, method = "bilinear") {
  raster::projectRaster(
    r_ll,
    crs = sp::CRS(crs_robinson),
    method = method
  )
}


# ----------------------------------------------------------
# Geometrías Natural Earth en Robinson
# ----------------------------------------------------------
get_land_robin <- function(r_proj) {
  
  land <- rnaturalearth::ne_download(
    scale = "medium",
    type = "land",
    category = "physical",
    returnclass = "sf"
  )
  
  land_rb <- sf::st_transform(land, crs = crs_robinson)
  
  sf::st_geometry(
    sf::st_crop(land_rb, sf::st_as_sfc(sf::st_bbox(r_proj)))
  )
}


get_coast_robin <- function(r_proj) {
  
  cst <- rnaturalearth::ne_download(
    scale = "medium",
    type = "coastline",
    category = "physical",
    returnclass = "sf"
  )
  
  cst_rb <- sf::st_transform(cst, crs = crs_robinson)
  
  sf::st_geometry(
    sf::st_crop(cst_rb, sf::st_as_sfc(sf::st_bbox(r_proj)))
  )
}


# ----------------------------------------------------------
# Paralelos y meridianos
# ----------------------------------------------------------
make_parallels_sf <- function(lats = c(-60, -40, 0, 40, 90),
                              step_deg = 0.5) {
  
  lines_list <- lapply(lats, function(la) {
    xx <- seq(-180, 180, by = step_deg)
    yy <- rep(la, length(xx))
    sf::st_linestring(cbind(xx, yy))
  })
  
  sf_obj <- sf::st_sfc(lines_list, crs = 4326)
  sf::st_transform(sf_obj, crs_robinson)
}


make_meridians_sf <- function(lons = c(-180, 180),
                              lat_min = -60,
                              lat_max = 90,
                              step_deg = 0.5,
                              r_proj = NULL) {
  
  lines_list <- lapply(lons, function(lo) {
    yy <- seq(lat_min, lat_max, by = step_deg)
    xx <- rep(lo, length(yy))
    sf::st_linestring(cbind(xx, yy))
  })
  
  sf_obj <- sf::st_sfc(lines_list, crs = 4326)
  sf_obj <- sf::st_transform(sf_obj, crs_robinson)
  
  if (!is.null(r_proj)) {
    sf_obj <- sf::st_crop(sf_obj, sf::st_as_sfc(sf::st_bbox(r_proj)))
  }
  
  sf_obj
}


# ----------------------------------------------------------
# Abrir dispositivo gráfico
# ----------------------------------------------------------
open_graphics <- function(outfile, width = 10, height = 5.8, bg = "white") {
  
  ext <- tolower(tools::file_ext(outfile))
  
  if (ext == "pdf") {
    
    grDevices::pdf(
      outfile,
      width = width,
      height = height,
      bg = bg,
      useDingbats = FALSE,
      onefile = FALSE
    )
    
  } else if (ext == "png") {
    
    grDevices::png(
      outfile,
      width = width,
      height = height,
      units = "in",
      res = 300,
      bg = bg
    )
    
  } else {
    
    warning("Extensión no reconocida. Se guardará como PDF.")
    
    grDevices::pdf(
      outfile,
      width = width,
      height = height,
      bg = bg,
      useDingbats = FALSE,
      onefile = FALSE
    )
  }
}


# ----------------------------------------------------------
# Eliminar ceros y valores negativos
# Para mapas de media anual
# ----------------------------------------------------------
positive_to_plot <- function(r, tol_zero = 0) {
  
  r[r <= tol_zero] <- NA
  r
}


# ----------------------------------------------------------
# Eliminar diferencias exactamente cero
# Para mapas de diferencias
# ----------------------------------------------------------
zero_to_NA <- function(r, tol_zero = 0) {
  
  if (tol_zero == 0) {
    r[r == 0] <- NA
  } else {
    r[abs(r) <= tol_zero] <- NA
  }
  
  r
}


# ==========================================================
# 3) Paletas discretas
# ==========================================================

# ----------------------------------------------------------
# Paleta de media anual de BA
# ----------------------------------------------------------
cols_annual <- c(
  "#346099", "#5984BD", "#88B0E3", "#A9E5E8",
  "#FFFF8C", "#FDAE61", "#D73027", "#A50026", "#870000"
)

breaks_annual <- c(
  0, 5, 10, 25, 50, 100, 250, 500, 750, Inf
)

labels_annual <- c(
  "< 5",
  "5-10",
  "10-25",
  "25-50",
  "50-100",
  "100-250",
  "250-500",
  "500-750",
  ">= 750"
)

palAnnual <- list(
  cols       = cols_annual,
  breaks     = breaks_annual,
  legend_at  = 1:9,
  legend_lab = labels_annual
)


# ----------------------------------------------------------
# Paleta de diferencias MRBA60 - producto
# ----------------------------------------------------------
cols_diff <- c(
  "#08306B", "#2171B5", "#6BAED6", "#A8DADC",
  "#FFFAC0",
  "#FDAEFB", "#DD65B6", "#C51B7D", "#7A0177"
)

breaks_diff <- c(
  -Inf, -75, -50, -25, -5, 5, 25, 50, 75, Inf
)

labels_diff <- c(
  "< -75",
  "-75 to -50",
  "-50 to -25",
  "-25 to -5",
  "-5 to 5",
  "5 to 25",
  "25 to 50",
  "50 to 75",
  ">= 75"
)

palDiff <- list(
  cols       = cols_diff,
  breaks     = breaks_diff,
  legend_at  = 1:9,
  legend_lab = labels_diff
)


# ==========================================================
# 4) Funciones de plotting
# ==========================================================

# ----------------------------------------------------------
# Mapa Robinson de media anual
# ----------------------------------------------------------
plot_robin_base <- function(r,
                            title,
                            subtitle_text,
                            pal,
                            out_path,
                            land,
                            coast,
                            parallels,
                            meridians,
                            legend_cex = 0.85) {
  
  open_graphics(out_path, width = 10, height = 5.8)
  par(mar = c(3.5, 0, 3, 0))
  
  # 0 y NA no se pintan
  r <- positive_to_plot(r, tol_zero = 0)
  
  # Clasificación:
  # [0,5), [5,10), [10,25), ..., [750, Inf)
  # Como r <= 0 ya es NA, el primer color representa 0 < x < 5.
  r_cut <- raster::calc(r, function(x) {
    cut(
      x,
      breaks = pal$breaks,
      labels = FALSE,
      include.lowest = FALSE,
      right = FALSE
    )
  })
  
  plot(
    land,
    col = "grey40",
    border = NA,
    main = title,
    cex.main = 2,
    line = 1,
    axes = FALSE
  )
  
  plot(
    r_cut,
    col = pal$cols,
    add = TRUE,
    legend = FALSE,
    zlim = c(1, length(pal$cols))
  )
  
  sub_expr <- bquote("Annual Mean BA " * km^2 ~ .(subtitle_text))
  mtext(sub_expr, side = 3, line = -0.6, cex = 1.2)
  
  plot(coast,     add = TRUE, border = "black", lwd = 0.8, lty = 1)
  plot(parallels, add = TRUE, col = "black",    lty = 1, lwd = 0.4)
  plot(meridians, add = TRUE, col = "black",    lty = 1, lwd = 0.4)
  
  fields::image.plot(
    zlim = c(1, length(pal$cols)),
    legend.only = TRUE,
    col = pal$cols,
    horizontal = TRUE,
    smallplot = c(0.15, 0.85, 0.08, 0.12),
    axis.args = list(
      at = pal$legend_at,
      labels = pal$legend_lab,
      cex.axis = legend_cex,
      lwd.ticks = 0,
      mgp = c(3, 0.2, 0)
    )
  )
  
  dev.off()
  
  message("Guardado: ", out_path)
}


# ----------------------------------------------------------
# Mapa Robinson de diferencias
# ----------------------------------------------------------
plot_robin_diff <- function(r,
                            title,
                            subtitle_text,
                            pal,
                            out_path,
                            land,
                            coast,
                            parallels,
                            meridians,
                            legend_cex = 0.75) {
  
  open_graphics(out_path, width = 10, height = 5.8)
  par(mar = c(3.5, 0, 3, 0))
  
  # Diferencia exactamente 0 no se pinta
  r <- zero_to_NA(r, tol_zero = 0)
  
  # Clasificación:
  # [-Inf,-75), [-75,-50), ..., [-5,5), [5,25), ..., [75,Inf)
  r_cut <- raster::calc(r, function(x) {
    cut(
      x,
      breaks = pal$breaks,
      labels = FALSE,
      include.lowest = TRUE,
      right = FALSE
    )
  })
  
  plot(
    land,
    col = "grey40",
    border = NA,
    main = title,
    cex.main = 2,
    line = 1,
    axes = FALSE
  )
  
  plot(
    r_cut,
    col = pal$cols,
    add = TRUE,
    legend = FALSE,
    zlim = c(1, length(pal$cols))
  )
  
  sub_expr <- bquote(italic(Delta) ~ "Annual Mean BA " * km^2 ~ .(subtitle_text))
  mtext(sub_expr, side = 3, line = -0.5, cex = 1.2)
  
  plot(coast,     add = TRUE, border = "black", lwd = 0.8, lty = 1)
  plot(parallels, add = TRUE, col = "black",    lty = 1, lwd = 0.4)
  plot(meridians, add = TRUE, col = "black",    lty = 1, lwd = 0.4)
  
  fields::image.plot(
    zlim = c(1, length(pal$cols)),
    legend.only = TRUE,
    col = pal$cols,
    horizontal = TRUE,
    smallplot = c(0.15, 0.85, 0.08, 0.12),
    axis.args = list(
      at = pal$legend_at,
      labels = pal$legend_lab,
      cex.axis = legend_cex,
      lwd.ticks = 0,
      mgp = c(3, 0.2, 0)
    )
  )
  
  dev.off()
  
  message("Guardado: ", out_path)
}


# ==========================================================
# 5) Cargar coordenadas
# ==========================================================
load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))

lon_vec <- if (!is.null(dim(lon))) {
  sort(unique(as.vector(lon)))
} else {
  as.vector(lon)
}

lat_vec <- if (!is.null(dim(lat))) {
  sort(unique(as.vector(lat)))
} else {
  as.vector(lat)
}

message("Longitudes: ", length(lon_vec))
message("Latitudes:  ", length(lat_vec))


# ==========================================================
# 6) Cargar datos de área quemada
# ==========================================================

# ----------------------------------------------------------
# MRBA60 armonizado
# Esperado: BA_MRBA60 en m2/mes
# Salida: BA_final en km2/mes
# ----------------------------------------------------------
ruta_RData_MRBA60 <- file.path(
  output_dir_RData,
  "MRBA60_BA_m2_monthly_2003_2024.RData"
)

load(ruta_RData_MRBA60)

if (!exists("BA_MRBA60")) {
  stop("No existe el objeto BA_MRBA60 dentro de MRBA60_BA_m2_monthly_2003_2024.RData.")
}

BA_final <- BA_MRBA60 / 1e6
BA_final[is.na(BA_final)] <- 0

message("Dimensiones MRBA60: ", paste(dim(BA_final), collapse = " x "))


# ----------------------------------------------------------
# FireCCI51
# Esperado: f51 en m2/mes
# Salida: BA_FireCCI51 en km2/mes
# ----------------------------------------------------------
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

if (!exists("f51")) {
  stop("No existe el objeto f51 dentro de FireCCI51_2003_2024_0.25degree.RData.")
}

BA_FireCCI51 <- f51 / 1e6
BA_FireCCI51[is.na(BA_FireCCI51)] <- 0

message("Dimensiones FireCCI51: ", paste(dim(BA_FireCCI51), collapse = " x "))


# ----------------------------------------------------------
# MCD64A1
# NetCDF mensual 2000-2024
# Se extrae 2003-2024: meses 27:290
# ----------------------------------------------------------
nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/MCD64CMQ_Monthly_2000-2024.nc"

nc_mcd <- ncdf4::nc_open(nc_path_mcd)
BA_MCD64 <- ncdf4::ncvar_get(nc_mcd, "band_data") / 1e6
ncdf4::nc_close(nc_mcd)

BA_MCD64 <- BA_MCD64[, , 27:290, drop = FALSE]
BA_MCD64[is.na(BA_MCD64)] <- 0

message("Dimensiones MCD64A1 2003-2024: ", paste(dim(BA_MCD64), collapse = " x "))


# ----------------------------------------------------------
# GFED5
# Esperado: burned_area en m2/mes o unidad equivalente según tu archivo
# Salida: BA_GFED5 en km2/mes
# ----------------------------------------------------------
nc_path_gfed <- "/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc"

nc_gfed <- ncdf4::nc_open(nc_path_gfed)
BA_GFED5 <- ncdf4::ncvar_get(nc_gfed, "burned_area") / 1e6
ncdf4::nc_close(nc_gfed)

BA_GFED5[is.na(BA_GFED5)] <- 0

message("Dimensiones GFED5: ", paste(dim(BA_GFED5), collapse = " x "))


# ==========================================================
# 7) Comprobaciones temporales
# ==========================================================

if (dim(BA_final)[3] < nmonths_0324) {
  stop("MRBA60 no tiene 264 meses.")
}

if (dim(BA_MCD64)[3] < nmonths_0324) {
  stop("MCD64A1 no tiene 264 meses para 2003-2024.")
}

if (dim(BA_GFED5)[3] < nmonths_0324) {
  stop("GFED5 no tiene 264 meses para 2003-2024.")
}

# FireCCI51 puede tener 240 o 264 meses según el archivo.
# El script detecta automáticamente su longitud.
nmonths_f51 <- dim(BA_FireCCI51)[3]

if (nmonths_f51 == 264) {
  
  dates_f51 <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
  text_f51  <- "(2003-2024)"
  message("FireCCI51 detectado como 2003-2024: 264 meses.")
  
} else if (nmonths_f51 == 240) {
  
  dates_f51 <- seq(as.Date("2003-01-01"), as.Date("2022-12-01"), by = "month")
  text_f51  <- "(2003-2022)"
  message("FireCCI51 detectado como 2003-2022: 240 meses.")
  
} else {
  
  stop(sprintf(
    "FireCCI51 tiene %d meses. Se esperaban 240 o 264.",
    nmonths_f51
  ))
}


# ==========================================================
# 8) Medias anuales por producto
# ==========================================================

# MRBA60, MCD64A1 y GFED5: siempre 2003-2024
BA60_mean_0324  <- monthly_to_annual_mean(BA_final, nmonths_0324)
MCD64_mean_0324 <- monthly_to_annual_mean(BA_MCD64, nmonths_0324)
GFED5_mean_0324 <- monthly_to_annual_mean(BA_GFED5, nmonths_0324)

# FireCCI51: periodo detectado automáticamente
BA51_mean <- monthly_to_annual_mean(BA_FireCCI51, nmonths_f51)

# MRBA60 recortado al periodo común con FireCCI51
BA60_mean_for51 <- monthly_to_annual_mean(BA_final, nmonths_f51)


# ==========================================================
# 9) Convertir medias a raster WGS84 y proyectar a Robinson
# ==========================================================

# ----------------------------------------------------------
# WGS84
# ----------------------------------------------------------
r60_ll   <- mat_to_raster_ll(BA60_mean_0324,  lon_vec, lat_vec)
r51_ll   <- mat_to_raster_ll(BA51_mean,       lon_vec, lat_vec)
rMCD_ll  <- mat_to_raster_ll(MCD64_mean_0324, lon_vec, lat_vec)
rGFED_ll <- mat_to_raster_ll(GFED5_mean_0324, lon_vec, lat_vec)


# ----------------------------------------------------------
# Robinson
# ----------------------------------------------------------
r60_rb   <- to_robin(r60_ll,   method = "bilinear")
r51_rb   <- to_robin(r51_ll,   method = "bilinear")
rMCD_rb  <- to_robin(rMCD_ll,  method = "bilinear")
rGFED_rb <- to_robin(rGFED_ll, method = "bilinear")


# ----------------------------------------------------------
# 0 o valores negativos no se pintan en mapas de media
# ----------------------------------------------------------
r60_rb   <- positive_to_plot(r60_rb,   tol_zero = 0)
r51_rb   <- positive_to_plot(r51_rb,   tol_zero = 0)
rMCD_rb  <- positive_to_plot(rMCD_rb,  tol_zero = 0)
rGFED_rb <- positive_to_plot(rGFED_rb, tol_zero = 0)


# ==========================================================
# 10) Geometrías en Robinson
# ==========================================================

land_geom      <- get_land_robin(r60_rb)
coast_geom     <- get_coast_robin(r60_rb)
parallels_geom <- make_parallels_sf(lats = c(-60, -40, 0, 40, 90))
meridians_geom <- make_meridians_sf(lons = c(-180, 180), r_proj = r60_rb)


# ==========================================================
# 11) Mapas de media anual
# ==========================================================

text_0324 <- "(2003-2024)"

plot_robin_base(
  r          = r60_rb,
  title      = "a) MRBA60",
  subtitle_text = text_0324,
  pal        = palAnnual,
  out_path   = file.path(output_dir_plot, "mapa_mean_BA_2003_2024_MRBA60_robin.pdf"),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)

plot_robin_base(
  r          = r51_rb,
  title      = "b) FireCCI51",
  subtitle_text = text_f51,
  pal        = palAnnual,
  out_path   = file.path(
    output_dir_plot,
    paste0("mapa_mean_BA_", format(min(dates_f51), "%Y"), "_", format(max(dates_f51), "%Y"), "_FireCCI51_robin.pdf")
  ),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)

plot_robin_base(
  r          = rMCD_rb,
  title      = "c) MCD64A1",
  subtitle_text = text_0324,
  pal        = palAnnual,
  out_path   = file.path(output_dir_plot, "mapa_mean_BA_2003_2024_MCD64A1_robin.pdf"),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)

plot_robin_base(
  r          = rGFED_rb,
  title      = "d) GFED5",
  subtitle_text = text_0324,
  pal        = palAnnual,
  out_path   = file.path(output_dir_plot, "mapa_mean_BA_2003_2024_GFED5_robin.pdf"),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)


# ==========================================================
# 12) Mapas de diferencias MRBA60 - producto
# ==========================================================

# ----------------------------------------------------------
# Diferencias de medias en el periodo común
# ----------------------------------------------------------

# MRBA60 - FireCCI51
# Periodo común detectado automáticamente: 2003-2022 o 2003-2024
diff_60_51 <- BA60_mean_for51 - BA51_mean

# MRBA60 - MCD64A1
# Periodo común: 2003-2024
diff_60_MCD <- BA60_mean_0324 - MCD64_mean_0324

# MRBA60 - GFED5
# Periodo común: 2003-2024
diff_60_GFED <- BA60_mean_0324 - GFED5_mean_0324


# ----------------------------------------------------------
# Diferencias a raster WGS84
# ----------------------------------------------------------
rD51_ll   <- mat_to_raster_ll(diff_60_51,   lon_vec, lat_vec)
rDMCD_ll  <- mat_to_raster_ll(diff_60_MCD,  lon_vec, lat_vec)
rDGFED_ll <- mat_to_raster_ll(diff_60_GFED, lon_vec, lat_vec)


# ----------------------------------------------------------
# Diferencias a Robinson
# ----------------------------------------------------------
rD51_rb   <- to_robin(rD51_ll,   method = "bilinear")
rDMCD_rb  <- to_robin(rDMCD_ll,  method = "bilinear")
rDGFED_rb <- to_robin(rDGFED_ll, method = "bilinear")


# ----------------------------------------------------------
# Diferencias exactamente 0 no se pintan
# Si quieres ocultar diferencias casi cero, usa tol_zero = 1e-12
# ----------------------------------------------------------
rD51_rb   <- zero_to_NA(rD51_rb,   tol_zero = 0)
rDMCD_rb  <- zero_to_NA(rDMCD_rb,  tol_zero = 0)
rDGFED_rb <- zero_to_NA(rDGFED_rb, tol_zero = 0)


# ----------------------------------------------------------
# Textos de periodo
# ----------------------------------------------------------
text_diff_51   <- text_f51
text_diff_0324 <- "(2003-2024)"


# ----------------------------------------------------------
# Plot diferencias
# ----------------------------------------------------------
plot_robin_diff(
  r          = rD51_rb,
  title      = "a) MRBA60 - FireCCI51",
  subtitle_text = text_diff_51,
  pal        = palDiff,
  out_path   = file.path(
    output_dir_plot,
    paste0(
      "mapa_diff_mean_BA_",
      format(min(dates_f51), "%Y"), "_",
      format(max(dates_f51), "%Y"),
      "_MRBA60_minus_FireCCI51_robin.pdf"
    )
  ),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)

plot_robin_diff(
  r          = rDMCD_rb,
  title      = "b) MRBA60 - MCD64A1",
  subtitle_text = text_diff_0324,
  pal        = palDiff,
  out_path   = file.path(output_dir_plot, "mapa_diff_mean_BA_2003_2024_MRBA60_minus_MCD64A1_robin.pdf"),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)

plot_robin_diff(
  r          = rDGFED_rb,
  title      = "c) MRBA60 - GFED5",
  subtitle_text = text_diff_0324,
  pal        = palDiff,
  out_path   = file.path(output_dir_plot, "mapa_diff_mean_BA_2003_2024_MRBA60_minus_GFED5_robin.pdf"),
  land       = land_geom,
  coast      = coast_geom,
  parallels  = parallels_geom,
  meridians  = meridians_geom
)


# ==========================================================
# 13) Guardar objetos principales usados en los mapas
# ==========================================================
# 
# save(
#   BA60_mean_0324,
#   BA51_mean,
#   MCD64_mean_0324,
#   GFED5_mean_0324,
#   diff_60_51,
#   diff_60_MCD,
#   diff_60_GFED,
#   dates_0324,
#   dates_f51,
#   palAnnual,
#   palDiff,
#   file = file.path(output_dir_RData, "Figure6_R1_mean_BA_and_differences_objects.RData")
# )

message("Proceso terminado correctamente.")
message("Figuras guardadas en: ", output_dir_plot)