# ==========================================================
# MAPA DE DIFERENCIA FIRECCIS311 - FIRECCI51
# Periodo: 2019-2024
#
# Diferencia:
#   FireCCIS311 - FireCCI51
#
# Unidad:
#   km2 yr-1 medio, 2019-2024
#
# Notas:
#   - La clase central -5 to 5 se pinta en amarillo claro.
#   - Solo las diferencias exactamente iguales a 0 se dejan sin color.
#   - Se excluye Antártida: lat > -60.
#   - Se añade borde exterior del mundo en proyección Robinson.
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
})


# ==========================================================
# Configuración de rutas
# ==========================================================
Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_base      <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_plot  <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS14_R1/"
output_dir_RData <- file.path(output_base, "RData")

dir.create(output_base,      showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_plot,  showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_RData, showWarnings = FALSE, recursive = TRUE)


# ==========================================================
# 1) Fechas
# ==========================================================
dates_0324 <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_1924 <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")

nmonths_0324 <- length(dates_0324)
nmonths_1924 <- length(dates_1924)

stopifnot(nmonths_0324 == 264)
stopifnot(nmonths_1924 == 72)


# ==========================================================
# 2) Funciones auxiliares
# ==========================================================

# ----------------------------------------------------------
# Detectar objeto cargado desde RData
# ----------------------------------------------------------
load_rdata_object <- function(path, possible_names) {
  
  if (!file.exists(path)) {
    stop("No existe el archivo: ", path)
  }
  
  env_tmp <- new.env()
  loaded_names <- load(path, envir = env_tmp)
  
  obj_name <- intersect(possible_names, loaded_names)
  
  if (length(obj_name) == 0) {
    stop(
      "No se encontró ninguno de los objetos esperados en:\n",
      path, "\n\n",
      "Objetos disponibles: ", paste(loaded_names, collapse = ", "), "\n\n",
      "Objetos esperados: ", paste(possible_names, collapse = ", ")
    )
  }
  
  get(obj_name[1], envir = env_tmp)
}


# ----------------------------------------------------------
# Conversión automática a km2
# Si los valores parecen estar en m2, divide entre 1e6.
# Si parecen estar ya en km2, los deja igual.
# ----------------------------------------------------------
to_km2_auto <- function(arr, product_name = "producto") {
  
  max_val <- suppressWarnings(max(arr, na.rm = TRUE))
  
  if (!is.finite(max_val)) {
    stop("No hay valores finitos en ", product_name)
  }
  
  if (max_val > 1e5) {
    message(product_name, ": valores detectados como m2 -> se convierten a km2.")
    arr <- arr / 1e6
  } else {
    message(product_name, ": valores detectados como km2 -> no se convierten.")
  }
  
  arr
}


# ----------------------------------------------------------
# Suma mensual -> suma anual -> media anual
# Devuelve matriz [lon, lat] en km2 yr-1
# ----------------------------------------------------------
monthly_to_annual_mean <- function(arr_3d) {
  
  if (length(dim(arr_3d)) != 3) {
    stop("Se esperaba un array 3D [lon, lat, time].")
  }
  
  nx <- dim(arr_3d)[1]
  ny <- dim(arr_3d)[2]
  nt <- dim(arr_3d)[3]
  
  if (nt %% 12 != 0) {
    stop("El número de meses debe ser múltiplo de 12.")
  }
  
  nyears <- nt / 12
  
  arr4 <- array(arr_3d, dim = c(nx, ny, 12, nyears))
  
  annual_sum  <- apply(arr4, c(1, 2, 4), sum, na.rm = TRUE)
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
# Geometría Natural Earth: land en Robinson
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


# ----------------------------------------------------------
# Geometría Natural Earth: coastline en Robinson
# ----------------------------------------------------------
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
# Paralelos
# ----------------------------------------------------------
make_parallels_sf <- function(lats = c(-60, -30, 0, 30, 60),
                              step_deg = 0.5) {
  
  lines_list <- lapply(lats, function(la) {
    xx <- seq(-180, 180, by = step_deg)
    yy <- rep(la, length(xx))
    sf::st_linestring(cbind(xx, yy))
  })
  
  sf_obj <- sf::st_sfc(lines_list, crs = 4326)
  sf::st_transform(sf_obj, crs_robinson)
}


# ----------------------------------------------------------
# Meridianos interiores
# ----------------------------------------------------------
make_meridians_sf <- function(lons = c(-120, -60, 0, 60, 120),
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
# Borde exterior del mapa Robinson
# ----------------------------------------------------------
make_robinson_frame_sf <- function(lat_min = -60,
                                   lat_max = 90,
                                   step_deg = 0.5) {
  
  yy_left <- seq(lat_min, lat_max, by = step_deg)
  xx_left <- rep(-180, length(yy_left))
  
  xx_top <- seq(-180, 180, by = step_deg)
  yy_top <- rep(lat_max, length(xx_top))
  
  yy_right <- seq(lat_max, lat_min, by = -step_deg)
  xx_right <- rep(180, length(yy_right))
  
  xx_bottom <- seq(180, -180, by = -step_deg)
  yy_bottom <- rep(lat_min, length(xx_bottom))
  
  coords <- rbind(
    cbind(xx_left, yy_left),
    cbind(xx_top, yy_top),
    cbind(xx_right, yy_right),
    cbind(xx_bottom, yy_bottom),
    cbind(-180, lat_min)
  )
  
  sf_obj <- sf::st_sfc(
    sf::st_linestring(coords),
    crs = 4326
  )
  
  sf::st_transform(sf_obj, crs_robinson)
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
    
  } else if (ext %in% c("jpg", "jpeg")) {
    
    grDevices::jpeg(
      outfile,
      width = width,
      height = height,
      units = "in",
      res = 300,
      quality = 100,
      bg = bg
    )
    
  } else {
    
    stop("Extensión no reconocida: ", ext)
  }
}


# ----------------------------------------------------------
# Eliminar diferencias exactamente cero
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
# 3) Paleta discreta de diferencias
# ==========================================================

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
# 4) Función de plot: diferencia en Robinson
# ==========================================================
plot_robin_diff <- function(r,
                            title,
                            subtitle_text,
                            pal,
                            out_path,
                            land,
                            coast,
                            parallels,
                            meridians,
                            frame,
                            legend_cex = 0.85) {
  
  open_graphics(out_path, width = 10, height = 5.8)
  par(mar = c(3.8, 0.2, 3.2, 0.2))
  
  # Solo diferencia exactamente 0 sin color.
  # Los valores entre -5 y 5 se pintan en amarillo claro.
  r <- zero_to_NA(r, tol_zero = 0)
  
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
    cex.main = 2.2,
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
  
  mtext(
    subtitle_text,
    side = 3,
    line = -0.6,
    cex = 1.35
  )
  
  plot(parallels, add = TRUE, col = "grey55", lty = 2, lwd = 0.5)
  plot(meridians, add = TRUE, col = "grey55", lty = 2, lwd = 0.5)
  plot(coast,     add = TRUE, border = "black", lwd = 0.75, lty = 1)
  plot(frame,     add = TRUE, col = "black", lwd = 1.1)
  
  fields::image.plot(
    zlim = c(1, length(pal$cols)),
    legend.only = TRUE,
    col = pal$cols,
    horizontal = TRUE,
    smallplot = c(0.13, 0.87, 0.075, 0.12),
    axis.args = list(
      at = pal$legend_at,
      labels = pal$legend_lab,
      cex.axis = legend_cex,
      lwd.ticks = 0,
      mgp = c(3, 0.25, 0)
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
# 6) Cargar FireCCI51
# ==========================================================
path_firecci51 <- file.path(
  dir_oss,
  "FireCCI51_2003_2024_0.25degree.RData"
)

BA_FireCCI51 <- load_rdata_object(
  path = path_firecci51,
  possible_names = c(
    "f51",
    "BA_FireCCI51",
    "BA_Fire51",
    "BA_Fire51_tot"
  )
)

BA_FireCCI51 <- to_km2_auto(BA_FireCCI51, "FireCCI51")
BA_FireCCI51[is.na(BA_FireCCI51)] <- 0

message("Dimensiones FireCCI51: ", paste(dim(BA_FireCCI51), collapse = " x "))

if (dim(BA_FireCCI51)[3] != 264) {
  stop(
    "FireCCI51 debe tener 264 meses para extraer 2019-2024. ",
    "Dimensión temporal encontrada: ", dim(BA_FireCCI51)[3]
  )
}

idx_1924_f51 <- which(dates_0324 %in% dates_1924)

if (length(idx_1924_f51) != 72) {
  stop("No se han encontrado correctamente los 72 meses de 2019-2024 en FireCCI51.")
}

BA_FireCCI51_1924 <- BA_FireCCI51[, , idx_1924_f51, drop = FALSE]

message("Dimensiones FireCCI51 2019-2024: ", paste(dim(BA_FireCCI51_1924), collapse = " x "))


# ==========================================================
# 7) Cargar FireCCIS311
# ==========================================================
path_fireccis311 <- file.path(
  dir_oss,
  "FireCCIS311_2019_2024_0.25degree.RData"
)

BA_FireCCIS311 <- load_rdata_object(
  path = path_fireccis311,
  possible_names = c(
    "s3",
    "BA_FireS3",
    "BA_FireS3_tot",
    "BA_FireCCIS311",
    "FireCCIS311"
  )
)

BA_FireCCIS311 <- to_km2_auto(BA_FireCCIS311, "FireCCIS311")
BA_FireCCIS311[is.na(BA_FireCCIS311)] <- 0

message("Dimensiones FireCCIS311: ", paste(dim(BA_FireCCIS311), collapse = " x "))

if (dim(BA_FireCCIS311)[3] == 72) {
  
  BA_FireCCIS311_1924 <- BA_FireCCIS311
  
} else if (dim(BA_FireCCIS311)[3] == 264) {
  
  BA_FireCCIS311_1924 <- BA_FireCCIS311[, , idx_1924_f51, drop = FALSE]
  
} else {
  
  stop(
    "FireCCIS311 debe tener 72 meses —2019-2024— o 264 meses —2003-2024—. ",
    "Dimensión temporal encontrada: ", dim(BA_FireCCIS311)[3]
  )
}

message("Dimensiones FireCCIS311 2019-2024: ", paste(dim(BA_FireCCIS311_1924), collapse = " x "))


# ==========================================================
# 8) Comprobar dimensiones espaciales
# ==========================================================
if (!all(dim(BA_FireCCI51_1924)[1:2] == dim(BA_FireCCIS311_1924)[1:2])) {
  stop(
    "Las dimensiones espaciales no coinciden:\n",
    "FireCCI51:   ", paste(dim(BA_FireCCI51_1924), collapse = " x "), "\n",
    "FireCCIS311: ", paste(dim(BA_FireCCIS311_1924), collapse = " x ")
  )
}

if (!all(dim(BA_FireCCI51_1924)[1:2] == c(length(lon_vec), length(lat_vec)))) {
  stop(
    "Las dimensiones espaciales no coinciden con lon/lat:\n",
    "Array: ", paste(dim(BA_FireCCI51_1924)[1:2], collapse = " x "), "\n",
    "lon/lat: ", length(lon_vec), " x ", length(lat_vec)
  )
}


# ==========================================================
# 9) Media anual 2019-2024
# ==========================================================
BA51_mean_1924 <- monthly_to_annual_mean(BA_FireCCI51_1924)
S311_mean_1924 <- monthly_to_annual_mean(BA_FireCCIS311_1924)


# ==========================================================
# 10) Diferencia FireCCIS311 - FireCCI51
# ==========================================================
diff_S311_51 <- S311_mean_1924 - BA51_mean_1924


# ==========================================================
# 11) Estadísticos rápidos de control
# ==========================================================
diff_values <- as.vector(diff_S311_51)
diff_values <- diff_values[is.finite(diff_values)]

summary_diff <- data.frame(
  Product_difference = "FireCCIS311 - FireCCI51",
  Period = "2019-2024",
  Unit = "km2 yr-1",
  n_cells = length(diff_values),
  mean_diff = mean(diff_values, na.rm = TRUE),
  median_diff = median(diff_values, na.rm = TRUE),
  sd_diff = sd(diff_values, na.rm = TRUE),
  min_diff = min(diff_values, na.rm = TRUE),
  p05_diff = as.numeric(quantile(diff_values, 0.05, na.rm = TRUE)),
  p25_diff = as.numeric(quantile(diff_values, 0.25, na.rm = TRUE)),
  p75_diff = as.numeric(quantile(diff_values, 0.75, na.rm = TRUE)),
  p95_diff = as.numeric(quantile(diff_values, 0.95, na.rm = TRUE)),
  max_diff = max(diff_values, na.rm = TRUE)
)

print(summary_diff)

write.csv(
  summary_diff,
  file = file.path(
    output_dir_plot,
    "summary_diff_mean_BA_2019_2024_FireCCIS311_minus_FireCCI51.csv"
  ),
  row.names = FALSE
)


# ==========================================================
# 12) Convertir diferencia a raster WGS84 y Robinson
# ==========================================================
rD_ll <- mat_to_raster_ll(
  mat = diff_S311_51,
  lon = lon_vec,
  lat = lat_vec,
  lat_min = -60
)

rD_rb <- to_robin(rD_ll, method = "bilinear")

rD_rb <- zero_to_NA(rD_rb, tol_zero = 0)


# ==========================================================
# 13) Geometrías en Robinson
# ==========================================================
land_geom <- get_land_robin(rD_rb)
coast_geom <- get_coast_robin(rD_rb)

parallels_geom <- make_parallels_sf(
  lats = c(-60, -30, 0, 30, 60)
)

meridians_geom <- make_meridians_sf(
  lons = c(-120, -60, 0, 60, 120),
  lat_min = -60,
  lat_max = 90,
  r_proj = rD_rb
)

frame_geom <- make_robinson_frame_sf(
  lat_min = -60,
  lat_max = 90
)


# ==========================================================
# 14) Guardar mapa en PDF, PNG y JPEG
# ==========================================================
title_map <- "FireCCIS311 - FireCCI51"

subtitle_map <- expression(
  Delta~"Annual Mean BA (km"^2~"yr"^-1*"), 2019-2024"
)

plot_robin_diff(
  r             = rD_rb,
  title         = title_map,
  subtitle_text = subtitle_map,
  pal           = palDiff,
  out_path      = file.path(
    output_dir_plot,
    "mapa_diff_mean_BA_2019_2024_FireCCIS311_minus_FireCCI51_robin.pdf"
  ),
  land          = land_geom,
  coast         = coast_geom,
  parallels     = parallels_geom,
  meridians     = meridians_geom,
  frame         = frame_geom,
  legend_cex    = 0.85
)

plot_robin_diff(
  r             = rD_rb,
  title         = title_map,
  subtitle_text = subtitle_map,
  pal           = palDiff,
  out_path      = file.path(
    output_dir_plot,
    "mapa_diff_mean_BA_2019_2024_FireCCIS311_minus_FireCCI51_robin.png"
  ),
  land          = land_geom,
  coast         = coast_geom,
  parallels     = parallels_geom,
  meridians     = meridians_geom,
  frame         = frame_geom,
  legend_cex    = 0.85
)

plot_robin_diff(
  r             = rD_rb,
  title         = title_map,
  subtitle_text = subtitle_map,
  pal           = palDiff,
  out_path      = file.path(
    output_dir_plot,
    "mapa_diff_mean_BA_2019_2024_FireCCIS311_minus_FireCCI51_robin.jpeg"
  ),
  land          = land_geom,
  coast         = coast_geom,
  parallels     = parallels_geom,
  meridians     = meridians_geom,
  frame         = frame_geom,
  legend_cex    = 0.85
)

# 
# # ==========================================================
# # 15) Guardar objetos principales
# # ==========================================================
# save(
#   BA51_mean_1924,
#   S311_mean_1924,
#   diff_S311_51,
#   summary_diff,
#   dates_1924,
#   palDiff,
#   file = file.path(
#     output_dir_RData,
#     "Diff_mean_BA_2019_2024_FireCCIS311_minus_FireCCI51.RData"
#   )
# )
# 
# 
# # ==========================================================
# # 16) Mensaje final
# # ==========================================================
# message("Proceso terminado correctamente.")
# message("Diferencia calculada como: FireCCIS311 - FireCCI51")
# message("Periodo: 2019-2024")
# message("Unidad: km2 yr-1")
# message("Figuras guardadas en: ", output_dir_plot)
# message("Objetos guardados en: ", output_dir_RData)