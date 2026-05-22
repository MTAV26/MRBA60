# ==========================================================
# BRASIL – ZOOM REGIONAL
# FireCCI51 / MRBA60H / MapBiomas / Diferencia
#
# - PLOTS individuales a–d + montaje 1×4
# - Inset overview Brasil + rectángulo del zoom
# - Inset en esquina inferior izquierda del panel a)
# - Zoom fijo latitudes: 15ºS a 5ºS
# - Paleta continua: rev(inferno(100))
# - Salida en:
#   /mnt/disco6tb/MRBA60/results/D1-Plots/Figure_S8_R1
#
# REQUIERE EN MEMORIA:
#   - lon
#   - lat
#   - BA_MRBA60H     [lon, lat, time] en km²/mes
#   - BA_FireCCI51   [lon, lat, time] en km²/mes
#   - template       SpatRaster global 0.25º
# ==========================================================

# rm(list = setdiff(ls(), c(
#   "lon", "lat",
#   "BA_MRBA60H", "BA_FireCCI51",
#   "template"
# )))

gc()

# ==========================================================
# 0. LIBRERÍAS
# ==========================================================

suppressPackageStartupMessages({
  library(terra)
  library(raster)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
  library(scales)
})

# ==========================================================
# 1. COMPROBACIONES INICIALES
# ==========================================================

required_objects <- c(
  "lon", "lat",
  "BA_MRBA60H", "BA_FireCCI51",
  "template"
)

missing_objects <- required_objects[!required_objects %in% ls()]

if (length(missing_objects) > 0) {
  stop(
    "Faltan objetos requeridos en memoria: ",
    paste(missing_objects, collapse = ", "),
    "\nCarga antes lon, lat, BA_MRBA60H, BA_FireCCI51 y template."
  )
}

cat("Objetos requeridos encontrados en memoria.\n")

cat("Dim BA_MRBA60H:", paste(dim(BA_MRBA60H), collapse = " x "), "\n")
cat("Dim BA_FireCCI51:", paste(dim(BA_FireCCI51), collapse = " x "), "\n")
cat("length(lon):", length(lon), "\n")
cat("length(lat):", length(lat), "\n")

stopifnot(dim(BA_MRBA60H)[1] == length(lon))
stopifnot(dim(BA_MRBA60H)[2] == length(lat))
stopifnot(dim(BA_FireCCI51)[1] == length(lon))
stopifnot(dim(BA_FireCCI51)[2] == length(lat))

# ==========================================================
# 2. PARÁMETROS / RUTAS
# ==========================================================

start_year_firecci <- 2003

ruta_MapBrasil <- "/mnt/disco6tb/Validacion/MapBrasil"

out_dir <- "/mnt/disco6tb/MRBA60/results/D1-Plots/Figure_S8_R1"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(out_dir)) {
  stop("No se ha podido crear la carpeta de salida: ", out_dir)
}

cat("Los archivos se guardarán en:\n", out_dir, "\n")

buffer_dist_m <- 25000

limit_diff_default <- 200

w_single_cm <- 12
h_single_cm <- 6
w_4pan_cm   <- 4 * w_single_cm
h_4pan_cm   <- h_single_cm

pal_ba <- rev(viridisLite::inferno(100))

pal_diff <- c(
  rev(viridisLite::inferno(50)),
  "#FFFFFF",
  viridisLite::inferno(50)
)

base_theme <- theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text  = element_text(color = "black"),
    axis.title = element_blank(),
    plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(0.35, "cm")
  )

# ==========================================================
# 3. REGIONES
# ==========================================================

regions <- list(
  list(
    name = "Amazonia_2005",
    year = 2005,
    lon  = c(-75, -55)
  )
)

# ==========================================================
# 4. PREPARAR LON/LAT
# ==========================================================

lon_glob <- lon
lat_glob <- lat

cat("Rango lon:", paste(range(lon_glob), collapse = " - "), "\n")
cat("Rango lat:", paste(range(lat_glob), collapse = " - "), "\n")

# ==========================================================
# 5. MÁSCARA BRASIL INTERIOR POR CENTROS DE CELDA
# ==========================================================

brasil_sf <- rnaturalearth::ne_countries(
  country = "Brazil",
  scale = "medium",
  returnclass = "sf"
)

brasil_sf <- sf::st_transform(brasil_sf, 4326)

brasil_proj <- sf::st_transform(brasil_sf, 3857)
brasil_union <- sf::st_union(brasil_proj)

brasil_inner_proj <- sf::st_buffer(
  brasil_union,
  dist = -buffer_dist_m
)

brasil_inner <- sf::st_transform(brasil_inner_proj, 4326)

# Recorte aproximado a Sudamérica para acelerar
ix_lon_sa <- which(lon_glob >= -90 & lon_glob <= -20)
ix_lat_sa <- which(lat_glob >= -40 & lat_glob <= 15)

lon_sa <- lon_glob[ix_lon_sa]
lat_sa <- lat_glob[ix_lat_sa]

grid_sa <- expand.grid(
  lon = lon_sa,
  lat = lat_sa
)

pts_sa <- sf::st_as_sf(
  grid_sa,
  coords = c("lon", "lat"),
  crs = 4326
)

inside_sa <- sf::st_within(
  pts_sa,
  brasil_inner,
  sparse = FALSE
)[, 1]

mask_sa <- matrix(
  inside_sa,
  nrow = length(lon_sa),
  ncol = length(lat_sa),
  byrow = FALSE
)

mask_brasil <- matrix(
  FALSE,
  nrow = length(lon_glob),
  ncol = length(lat_glob)
)

mask_brasil[ix_lon_sa, ix_lat_sa] <- mask_sa

cat("Celdas dentro de Brasil interior:", sum(mask_brasil, na.rm = TRUE), "\n")

# ==========================================================
# 6. CARGA MAPBIOMAS -> ARRAY [lon, lat, time]
# ==========================================================

patron_mb <- "^burned_km2_025deg_(\\d{4})\\.tif$"

archivos_mb <- list.files(
  path = ruta_MapBrasil,
  pattern = patron_mb,
  full.names = TRUE
)

if (length(archivos_mb) == 0) {
  stop("No se encontraron GeoTIFF MapBiomas en: ", ruta_MapBrasil)
}

archivos_mb <- archivos_mb[order(archivos_mb)]

years_mb <- as.numeric(
  sub(patron_mb, "\\1", basename(archivos_mb))
)

cat("Años MapBiomas disponibles:\n")
print(years_mb)

stk_mb <- raster::stack(archivos_mb)

e_mb  <- raster::extent(stk_mb)
dim_x <- raster::ncol(stk_mb)
dim_y <- raster::nrow(stk_mb)
nt_mb <- raster::nlayers(stk_mb)

ix_lon <- which(lon_glob >= e_mb@xmin & lon_glob <= e_mb@xmax)
ix_lat <- which(lat_glob >= e_mb@ymin & lat_glob <= e_mb@ymax)

lon_br <- lon_glob[ix_lon]
lat_br <- lat_glob[ix_lat]

mask_br_subset <- mask_brasil[ix_lon, ix_lat]

cat("Dim MapBiomas raster:", dim_x, "x", dim_y, "x", nt_mb, "\n")
cat("Dim lon_br/lat_br:", length(lon_br), "x", length(lat_br), "\n")

BA_MapBrasil <- array(
  NA_real_,
  dim = c(length(lon_br), length(lat_br), nt_mb)
)

for (i in seq_len(nt_mb)) {
  
  cat("Leyendo MapBiomas capa", i, "año", years_mb[i], "\n")
  
  mat <- raster::as.matrix(stk_mb[[i]], na.rm = FALSE)
  
  if (!is.matrix(mat)) {
    mat <- matrix(
      mat,
      nrow = dim_y,
      ncol = dim_x,
      byrow = TRUE
    )
  } else if (nrow(mat) != dim_y || ncol(mat) != dim_x) {
    mat <- matrix(
      as.vector(mat),
      nrow = dim_y,
      ncol = dim_x,
      byrow = TRUE
    )
  }
  
  dimnames(mat) <- NULL
  
  # raster::as.matrix devuelve filas norte -> sur.
  # Se pasa a matriz [lon, lat].
  mat_flip <- mat[dim_y:1, , drop = FALSE]
  slice <- t(mat_flip)
  dimnames(slice) <- NULL
  
  if (nrow(slice) != length(lon_br) || ncol(slice) != length(lat_br)) {
    stop(
      "Dimensiones MapBiomas no cuadran con lon_br/lat_br. ",
      "Revisa rejilla/extent."
    )
  }
  
  slice[!mask_br_subset] <- NA
  BA_MapBrasil[, , i] <- slice
}

# ==========================================================
# 7. HELPERS
# ==========================================================

get_mapbiomas_annual_mat <- function(year_target) {
  
  if (!(year_target %in% years_mb)) {
    stop(
      "Año ", year_target,
      " no está en MapBiomas. Disponibles: ",
      paste(years_mb, collapse = ", ")
    )
  }
  
  BA_MapBrasil[, , which(years_mb == year_target)]
}

get_firecci_annual_mat <- function(year_target, BA_FireCCI) {
  
  start_idx <- (year_target - start_year_firecci) * 12 + 1
  idx <- start_idx:(start_idx + 11)
  
  if (min(idx) < 1 || max(idx) > dim(BA_FireCCI)[3]) {
    stop("Índices temporales fuera del array para el año ", year_target)
  }
  
  arr_year_glob <- BA_FireCCI[, , idx, drop = FALSE]
  arr_year_br   <- arr_year_glob[ix_lon, ix_lat, , drop = FALSE]
  
  mat_annual <- apply(
    arr_year_br,
    c(1, 2),
    sum,
    na.rm = TRUE
  )
  
  mat_annual[!mask_br_subset] <- NA
  
  mat_annual
}

mat_lonlat_to_rast <- function(mat_lonlat, lon_vec, lat_vec, template_crs) {
  
  if (!is.matrix(mat_lonlat)) {
    stop("mat_lonlat debe ser una matriz.")
  }
  
  if (nrow(mat_lonlat) != length(lon_vec)) {
    stop("nrow(mat_lonlat) no coincide con length(lon_vec).")
  }
  
  if (ncol(mat_lonlat) != length(lat_vec)) {
    stop("ncol(mat_lonlat) no coincide con length(lat_vec).")
  }
  
  mat_latlon <- t(mat_lonlat)
  
  # terra espera fila 1 = norte.
  mat_flip <- mat_latlon[nrow(mat_latlon):1, , drop = FALSE]
  
  d_lon <- abs(lon_vec[2] - lon_vec[1])
  d_lat <- abs(lat_vec[2] - lat_vec[1])
  
  ext_br <- terra::ext(
    min(lon_vec) - d_lon / 2,
    max(lon_vec) + d_lon / 2,
    min(lat_vec) - d_lat / 2,
    max(lat_vec) + d_lat / 2
  )
  
  r <- terra::rast(
    ncols = length(lon_vec),
    nrows = length(lat_vec),
    extent = ext_br,
    crs = template_crs
  )
  
  terra::values(r) <- as.vector(t(mat_flip))
  
  r
}

make_overview_inset <- function(br_sf, xlim_zoom, ylim_zoom) {
  
  ggplot() +
    geom_sf(
      data = br_sf,
      fill = "grey90",
      color = "grey0",
      linewidth = 0.3
    ) +
    geom_rect(
      aes(
        xmin = xlim_zoom[1],
        xmax = xlim_zoom[2],
        ymin = ylim_zoom[1],
        ymax = ylim_zoom[2]
      ),
      fill = NA,
      color = "red",
      linewidth = 0.6
    ) +
    coord_sf(
      crs = sf::st_crs(4326),
      expand = FALSE
    ) +
    theme_void()
}

# ==========================================================
# 8. FUNCIÓN PRINCIPAL DE PLOT Y GUARDADO
# ==========================================================

plot_zoom_and_save <- function(region,
                               F51_r,
                               F60_r,
                               MB_r,
                               DIFF_r,
                               br_sf,
                               out_dir,
                               limit_diff = limit_diff_default,
                               w_single = w_single_cm,
                               h_single = h_single_cm,
                               w_4pan = w_4pan_cm,
                               h_4pan = h_4pan_cm) {
  
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  if (!dir.exists(out_dir)) {
    stop("No existe la carpeta de salida dentro de plot_zoom_and_save(): ", out_dir)
  }
  
  xlim_zoom <- region$lon
  ylim_zoom <- c(-15, -5)
  
  ext_zoom <- terra::ext(
    xlim_zoom[1],
    xlim_zoom[2],
    ylim_zoom[1],
    ylim_zoom[2]
  )
  
  F51z <- terra::crop(F51_r, ext_zoom)
  F60z <- terra::crop(F60_r, ext_zoom)
  MBz  <- terra::crop(MB_r,  ext_zoom)
  DIFz <- terra::crop(DIFF_r, ext_zoom)
  
  F51z[F51z == 0] <- NA
  F60z[F60z == 0] <- NA
  MBz[MBz == 0]   <- NA
  DIFz[DIFz == 0] <- NA
  
  df_51 <- as.data.frame(F51z, xy = TRUE, na.rm = TRUE)
  df_60 <- as.data.frame(F60z, xy = TRUE, na.rm = TRUE)
  df_mb <- as.data.frame(MBz,  xy = TRUE, na.rm = TRUE)
  df_df <- as.data.frame(DIFz, xy = TRUE, na.rm = TRUE)
  
  if (ncol(df_51) < 3) stop("df_51 no tiene datos para el zoom.")
  if (ncol(df_60) < 3) stop("df_60 no tiene datos para el zoom.")
  if (ncol(df_mb) < 3) stop("df_mb no tiene datos para el zoom.")
  if (ncol(df_df) < 3) stop("df_df no tiene datos para el zoom.")
  
  names(df_51)[3] <- "BA"
  names(df_60)[3] <- "BA"
  names(df_mb)[3] <- "BA"
  names(df_df)[3] <- "diff_BA"
  
  all_ba <- c(df_51$BA, df_60$BA, df_mb$BA)
  all_ba <- all_ba[is.finite(all_ba)]
  
  cap_ba <- as.numeric(
    stats::quantile(
      all_ba,
      probs = 0.995,
      na.rm = TRUE
    )
  )
  
  if (!is.finite(cap_ba) || cap_ba <= 0) {
    cap_ba <- max(all_ba, na.rm = TRUE)
  }
  
  if (is.na(limit_diff)) {
    rng <- range(df_df$diff_BA, na.rm = TRUE)
    limit_diff <- max(abs(rng))
  }
  
  br_sf_zoom <- sf::st_crop(
    br_sf,
    xmin = xlim_zoom[1],
    xmax = xlim_zoom[2],
    ymin = ylim_zoom[1],
    ymax = ylim_zoom[2]
  )
  
  lon_ticks <- seq(
    floor(xlim_zoom[1] / 5) * 5,
    ceiling(xlim_zoom[2] / 5) * 5,
    by = 5
  )
  
  lat_ticks <- seq(
    floor(ylim_zoom[1] / 2) * 2,
    ceiling(ylim_zoom[2] / 2) * 2,
    by = 2
  )
  
  title_suffix <- paste0(region$name, " – ", region$year)
  
  p_over <- make_overview_inset(
    br_sf = br_sf,
    xlim_zoom = xlim_zoom,
    ylim_zoom = ylim_zoom
  )
  
  inset_pos <- list(
    left = 0.02,
    bottom = 0.02,
    right = 0.30,
    top = 0.35
  )
  
  add_inset <- function(p) {
    p + patchwork::inset_element(
      p_over,
      left   = inset_pos$left,
      bottom = inset_pos$bottom,
      right  = inset_pos$right,
      top    = inset_pos$top
    )
  }
  
  p_a <- ggplot() +
    geom_sf(
      data = br_sf_zoom,
      fill = "grey80",
      color = "grey0",
      linewidth = 0.3
    ) +
    geom_tile(
      data = df_51,
      aes(x = x, y = y, fill = BA)
    ) +
    scale_fill_gradientn(
      colours = pal_ba,
      limits = c(0, cap_ba),
      oob = scales::squish,
      name = paste0("Annual BA ", region$year, " (km²)")
    ) +
    scale_x_continuous("", breaks = lon_ticks, expand = c(0, 0)) +
    scale_y_continuous("", breaks = lat_ticks, expand = c(0, 0)) +
    coord_sf(
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      expand = FALSE,
      crs = sf::st_crs(4326)
    ) +
    base_theme +
    ggtitle(paste0("a) FireCCI51 (", title_suffix, ")"))
  
  p_b <- ggplot() +
    geom_sf(
      data = br_sf_zoom,
      fill = "grey80",
      color = "grey0",
      linewidth = 0.3
    ) +
    geom_tile(
      data = df_60,
      aes(x = x, y = y, fill = BA)
    ) +
    scale_fill_gradientn(
      colours = pal_ba,
      limits = c(0, cap_ba),
      oob = scales::squish,
      name = paste0("Annual BA ", region$year, " (km²)")
    ) +
    scale_x_continuous("", breaks = lon_ticks, expand = c(0, 0)) +
    scale_y_continuous("", breaks = lat_ticks, expand = c(0, 0)) +
    coord_sf(
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      expand = FALSE,
      crs = sf::st_crs(4326)
    ) +
    base_theme +
    ggtitle(paste0("b) MRBA60H (", title_suffix, ")"))
  
  p_c <- ggplot() +
    geom_sf(
      data = br_sf_zoom,
      fill = "grey80",
      color = "grey0",
      linewidth = 0.3
    ) +
    geom_tile(
      data = df_mb,
      aes(x = x, y = y, fill = BA)
    ) +
    scale_fill_gradientn(
      colours = pal_ba,
      limits = c(0, cap_ba),
      oob = scales::squish,
      name = paste0("Annual BA ", region$year, " (km²)")
    ) +
    scale_x_continuous("", breaks = lon_ticks, expand = c(0, 0)) +
    scale_y_continuous("", breaks = lat_ticks, expand = c(0, 0)) +
    coord_sf(
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      expand = FALSE,
      crs = sf::st_crs(4326)
    ) +
    base_theme +
    ggtitle(paste0("c) MapBiomas (", title_suffix, ")"))
  
  p_d <- ggplot() +
    geom_sf(
      data = br_sf_zoom,
      fill = "grey80",
      color = "grey0",
      linewidth = 0.3
    ) +
    geom_tile(
      data = df_df,
      aes(x = x, y = y, fill = diff_BA)
    ) +
    scale_fill_gradientn(
      colours = pal_diff,
      limits = c(-limit_diff, limit_diff),
      oob = scales::squish,
      name = expression(Delta~"BA (km"^2*")"~"MRBA60H - MapBiomas")
    ) +
    scale_x_continuous("", breaks = lon_ticks, expand = c(0, 0)) +
    scale_y_continuous("", breaks = lat_ticks, expand = c(0, 0)) +
    coord_sf(
      xlim = xlim_zoom,
      ylim = ylim_zoom,
      expand = FALSE,
      crs = sf::st_crs(4326)
    ) +
    base_theme +
    ggtitle(paste0("d) Difference (", title_suffix, ")"))
  
  p_a <- add_inset(p_a)
  
  file_a <- file.path(
    out_dir,
    paste0("Zoom_", region$name, "_a_F51_", region$year, "_cont.pdf")
  )
  
  file_b <- file.path(
    out_dir,
    paste0("Zoom_", region$name, "_b_F60_", region$year, "_cont.pdf")
  )
  
  file_c <- file.path(
    out_dir,
    paste0("Zoom_", region$name, "_c_MB_", region$year, "_cont.pdf")
  )
  
  file_d <- file.path(
    out_dir,
    paste0("Zoom_", region$name, "_d_Diff_", region$year, "_cont.pdf")
  )
  
  ggsave(
    filename = file_a,
    plot = p_a,
    width = w_single,
    height = h_single,
    units = "cm"
  )
  
  ggsave(
    filename = file_b,
    plot = p_b,
    width = w_single,
    height = h_single,
    units = "cm"
  )
  
  ggsave(
    filename = file_c,
    plot = p_c,
    width = w_single,
    height = h_single,
    units = "cm"
  )
  
  ggsave(
    filename = file_d,
    plot = p_d,
    width = w_single,
    height = h_single,
    units = "cm"
  )
  
  p_all <- (p_a + p_b + p_c + p_d) +
    patchwork::plot_layout(ncol = 4, nrow = 1)
  
  file_all <- file.path(
    out_dir,
    paste0("Zoom_", region$name, "_4panels_", region$year, "_cont.pdf")
  )
  
  ggsave(
    filename = file_all,
    plot = p_all,
    width = w_4pan,
    height = h_4pan,
    units = "cm"
  )
  
  files_out <- c(
    a = file_a,
    b = file_b,
    c = file_c,
    d = file_d,
    all = file_all
  )
  
  cat("\nArchivos guardados:\n")
  print(files_out)
  
  cat("\nComprobación file.exists():\n")
  print(file.exists(files_out))
  
  list(
    p_a = p_a,
    p_b = p_b,
    p_c = p_c,
    p_d = p_d,
    p_all = p_all,
    cap_ba = cap_ba,
    limit_diff = limit_diff,
    files = files_out
  )
}

# ==========================================================
# 9. LOOP REGIONES
# ==========================================================

br_sf <- rnaturalearth::ne_countries(
  country = "Brazil",
  scale = "medium",
  returnclass = "sf"
)

br_sf <- sf::st_transform(br_sf, 4326)

results_zoom <- list()

for (region in regions) {
  
  year_target <- region$year
  
  cat("\n==============================\n")
  cat("Procesando:", region$name, "(", year_target, ")\n")
  cat("==============================\n")
  
  MB_mat <- get_mapbiomas_annual_mat(year_target)
  MB_rast <- mat_lonlat_to_rast(
    mat_lonlat = MB_mat,
    lon_vec = lon_br,
    lat_vec = lat_br,
    template_crs = terra::crs(template)
  )
  
  F60_mat <- get_firecci_annual_mat(
    year_target = year_target,
    BA_FireCCI = BA_MRBA60H
  )
  
  F51_mat <- get_firecci_annual_mat(
    year_target = year_target,
    BA_FireCCI = BA_FireCCI51
  )
  
  F60_rast <- mat_lonlat_to_rast(
    mat_lonlat = F60_mat,
    lon_vec = lon_br,
    lat_vec = lat_br,
    template_crs = terra::crs(template)
  )
  
  F51_rast <- mat_lonlat_to_rast(
    mat_lonlat = F51_mat,
    lon_vec = lon_br,
    lat_vec = lat_br,
    template_crs = terra::crs(template)
  )
  
  DIFF_rast <- F60_rast - MB_rast
  DIFF_rast[DIFF_rast == 0] <- NA
  
  cat("\n--- SUMA TOTAL BRASIL INTERIOR", year_target, "(km²) ---\n")
  cat("MapBiomas :", sum(MB_mat, na.rm = TRUE), "\n")
  cat("MRBA60H   :", sum(F60_mat, na.rm = TRUE), "\n")
  cat("FireCCI51 :", sum(F51_mat, na.rm = TRUE), "\n")
  
  res <- plot_zoom_and_save(
    region = region,
    F51_r = F51_rast,
    F60_r = F60_rast,
    MB_r = MB_rast,
    DIFF_r = DIFF_rast,
    br_sf = br_sf,
    out_dir = out_dir,
    limit_diff = limit_diff_default
  )
  
  results_zoom[[region$name]] <- res
  
  print(res$p_all)
  
  cat("\ncap_ba p99.5 usado:", res$cap_ba, "\n")
  cat("limit_diff usado:", res$limit_diff, "\n")
  cat("Guardados en:\n")
  print(res$files)
}

cat("\nListo.\n")
cat("Carpeta final de salida:\n", out_dir, "\n")
