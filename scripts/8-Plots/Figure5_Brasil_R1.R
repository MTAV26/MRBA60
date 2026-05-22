rm(list = ls())
graphics.off()
gc()

# ==========================================================
# 0. LIBRERÍAS
# ==========================================================
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggplot2)
library(sf)
library(dplyr)
library(scales)
library(patchwork)

sf::sf_use_s2(FALSE)

# ==========================================================
# 1. PARÁMETROS
# ==========================================================
year_target <- 2005
start_year  <- 2003

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"
output_base <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_RData <- file.path(output_base, "RData")

path_mapbiomas <- "/mnt/disco6tb/Validacion/MapBrasil"

dir_fig <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure5_R1"
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

# Distancia para evitar celdas de borde de Brasil
buffer_dist_m <- 25000

# Proyección métrica recomendada para Brasil
# SIRGAS 2000 / Brazil Polyconic
crs_brazil_m <- 5880

# ==========================================================
# 2. CARGAR MAPBIOMAS BRASIL 2005
#    Archivos esperados: burned_km2_025deg_YYYY.tif
# ==========================================================
pattern_mb <- "^burned_km2_025deg_(\\d{4})\\.tif$"

files_mb <- sort(list.files(
  path       = path_mapbiomas,
  pattern    = pattern_mb,
  full.names = TRUE
))

if (length(files_mb) == 0) {
  stop("No se encontraron archivos MapBiomas en: ", path_mapbiomas)
}

years_mb <- as.numeric(sub(pattern_mb, "\\1", basename(files_mb)))

cat("Años disponibles en MapBiomas:\n")
print(years_mb)

if (!(year_target %in% years_mb)) {
  stop("El año ", year_target, " no está disponible en MapBiomas.")
}

file_mb_2005 <- files_mb[years_mb == year_target]

cat("Archivo MapBiomas seleccionado:\n")
print(file_mb_2005)

MB_Brazil_annual_raw <- rast(file_mb_2005)
names(MB_Brazil_annual_raw) <- paste0("MapBiomas_", year_target)

cat("\nMapBiomas Brasil ", year_target, "\n")
print(MB_Brazil_annual_raw)

cat("Suma MapBiomas Brasil completa antes de recorte (km²):",
    global(MB_Brazil_annual_raw, "sum", na.rm = TRUE)$sum, "\n")

# ==========================================================
# 3. CARGAR LAT/LON Y ARRAYS MRBA60 / FireCCI51
# ==========================================================
load(file.path(dir_oss, "longitude.RData"))
load(file.path(dir_oss, "latitude.RData"))

cat("\nRango lat:", paste(range(lat), collapse = " - "), "\n")
cat("Rango lon:", paste(range(lon), collapse = " - "), "\n")

# ----------------------------------------------------------
# MRBA60 armonizado
# ----------------------------------------------------------
ruta_RData <- file.path(
  output_dir_RData,
  "MRBA60_BA_m2_monthly_2003_2024.RData"
)

load(ruta_RData)

# BA_MRBA60 está en m². Convertimos a km².
BA_MRBA60 <- BA_MRBA60 / 1e6
BA_MRBA60[is.na(BA_MRBA60)] <- 0

# ----------------------------------------------------------
# FireCCI51
# ----------------------------------------------------------
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

BA_FireCCI51 <- f51 / 1e6
BA_FireCCI51[BA_FireCCI51 == 0] <- NA

rm(f51)
gc()

cat("\nDims MRBA60:", paste(dim(BA_MRBA60), collapse = " x "), "\n")
cat("Dims FireCCI51:", paste(dim(BA_FireCCI51), collapse = " x "), "\n")

cat("length(lon):", length(lon), " | dim X arrays:", dim(BA_MRBA60)[1], "\n")
cat("length(lat):", length(lat), " | dim Y arrays:", dim(BA_MRBA60)[2], "\n")

stopifnot(dim(BA_MRBA60)[1] == length(lon))
stopifnot(dim(BA_MRBA60)[2] == length(lat))
stopifnot(all(dim(BA_MRBA60) == dim(BA_FireCCI51)))

# ==========================================================
# 4. ÍNDICES DE LOS 12 MESES DEL AÑO OBJETIVO
#    Serie mensual desde 2003-01
# ==========================================================
start_idx <- (year_target - start_year) * 12 + 1
end_idx   <- start_idx + 11
idx_year  <- start_idx:end_idx

cat("\nÍndices", year_target, ":", paste(idx_year, collapse = ","), "\n")

arr60_year <- BA_MRBA60[, , idx_year, drop = FALSE]
arr51_year <- BA_FireCCI51[, , idx_year, drop = FALSE]

cat("Dim arr60_year:", paste(dim(arr60_year), collapse = " x "), "\n")
cat("Dim arr51_year:", paste(dim(arr51_year), collapse = " x "), "\n")

cat("Rango arr60_year:", paste(range(arr60_year, na.rm = TRUE), collapse = " - "), "\n")
cat("Rango arr51_year:", paste(range(arr51_year, na.rm = TRUE), collapse = " - "), "\n")

cat("Suma global array MRBA60", year_target, "(km²):",
    sum(arr60_year, na.rm = TRUE), "\n")

cat("Suma global array FireCCI51", year_target, "(km²):",
    sum(arr51_year, na.rm = TRUE), "\n")

# ==========================================================
# 5. CREAR TEMPLATE GLOBAL 0.25º DESDE LON/LAT
# ==========================================================
d_lon <- abs(lon[2] - lon[1])
d_lat <- abs(lat[2] - lat[1])

ext_glob <- ext(
  min(lon) - d_lon / 2,
  max(lon) + d_lon / 2,
  min(lat) - d_lat / 2,
  max(lat) + d_lat / 2
)

template <- rast(
  ncols  = length(lon),
  nrows  = length(lat),
  extent = ext_glob,
  crs    = "EPSG:4326"
)

ncols <- ncol(template)
nrows <- nrow(template)

cat("\ntemplate ncols/nrows:", ncols, nrows, "\n")
cat("length(lon)/length(lat):", length(lon), length(lat), "\n")

stopifnot(ncols == length(lon))
stopifnot(nrows == length(lat))

# ==========================================================
# 6. RASTERIZAR MRBA60 2005 CON GEOMETRÍA GLOBAL
# ==========================================================
stack60_year <- rast(
  ncols  = ncols,
  nrows  = nrows,
  nlyrs  = 12,
  extent = ext_glob,
  crs    = crs(template)
)

for (m in 1:12) {
  cat("Rellenando MRBA60 mes", m, "\n")
  
  slice_2d <- arr60_year[, , m]
  mat_latlon <- t(slice_2d)
  mat_latlon_flip <- mat_latlon[nrow(mat_latlon):1, ]
  
  stack60_year[[m]][] <- as.vector(t(mat_latlon_flip))
}

names(stack60_year) <- sprintf("MRBA60_%s_%02d", year_target, 1:12)

cat("Suma global raster MRBA60", year_target, "(km²):",
    global(sum(stack60_year, na.rm = TRUE), "sum", na.rm = TRUE)$sum, "\n")

# ==========================================================
# 7. RASTERIZAR FireCCI51 2005 CON LA MISMA GEOMETRÍA
# ==========================================================
stack51_year <- rast(
  ncols  = ncols,
  nrows  = nrows,
  nlyrs  = 12,
  extent = ext_glob,
  crs    = crs(template)
)

for (m in 1:12) {
  cat("Rellenando FireCCI51 mes", m, "\n")
  
  slice_2d <- arr51_year[, , m]
  mat_latlon <- t(slice_2d)
  mat_latlon_flip <- mat_latlon[nrow(mat_latlon):1, ]
  
  stack51_year[[m]][] <- as.vector(t(mat_latlon_flip))
}

names(stack51_year) <- sprintf("FireCCI51_%s_%02d", year_target, 1:12)

cat("Suma global raster FireCCI51", year_target, "(km²):",
    global(sum(stack51_year, na.rm = TRUE), "sum", na.rm = TRUE)$sum, "\n")

# ==========================================================
# 8. PREPARAR BRASIL INTERIOR CON BUFFER MÉTRICO REAL
# ==========================================================
brazil_sf_full <- ne_countries(
  country     = "Brazil",
  scale       = "large",
  returnclass = "sf"
)

brazil_sf_full <- brazil_sf_full |>
  st_transform(4326) |>
  st_make_valid()

brazil_proj <- brazil_sf_full |>
  st_transform(crs_brazil_m) |>
  st_make_valid()

brazil_union <- brazil_proj |>
  st_union() |>
  st_make_valid()

brazil_inner_proj <- st_buffer(
  brazil_union,
  dist = -buffer_dist_m
) |>
  st_make_valid()

if (st_is_empty(brazil_inner_proj)) {
  stop("El buffer negativo ha generado una geometría vacía. Revisa buffer_dist_m.")
}

brazil_inner_sf <- brazil_inner_proj |>
  st_transform(4326) |>
  st_make_valid()

brazil_inner_vect <- terra::vect(brazil_inner_sf)

brazil_full_plot  <- brazil_sf_full
brazil_inner_plot <- st_as_sf(brazil_inner_vect)

bbox_br <- st_bbox(brazil_sf_full)

xlim_br <- c(bbox_br["xmin"], bbox_br["xmax"])
ylim_br <- c(bbox_br["ymin"], bbox_br["ymax"])

# Longitud cada 10 grados
lon_ticks <- seq(
  floor(xlim_br[1] / 10) * 10,
  ceiling(xlim_br[2] / 10) * 10,
  by = 10
)

# Latitud cada 5 grados
lat_ticks <- seq(
  floor(ylim_br[1] / 5) * 5,
  ceiling(ylim_br[2] / 5) * 5,
  by = 5
)

# ----------------------------------------------------------
# 8.1 CHECK NUMÉRICO DEL BUFFER
# ----------------------------------------------------------
area_brazil_full_km2  <- as.numeric(st_area(brazil_union)) / 1e6
area_brazil_inner_km2 <- as.numeric(st_area(brazil_inner_proj)) / 1e6
area_removed_km2      <- area_brazil_full_km2 - area_brazil_inner_km2

cat("\n--- CHECK BUFFER BRASIL ---\n")
cat("CRS métrico usado:", crs_brazil_m, "\n")
cat("Buffer negativo aplicado (m):", buffer_dist_m, "\n")
cat("Área Brasil completo km²:", area_brazil_full_km2, "\n")
cat("Área Brasil interior -25 km km²:", area_brazil_inner_km2, "\n")
cat("Área eliminada por buffer km²:", area_removed_km2, "\n")
cat("Porcentaje eliminado (%):", 100 * area_removed_km2 / area_brazil_full_km2, "\n")

# ----------------------------------------------------------
# 8.2 CHECK VISUAL DEL BUFFER
# ----------------------------------------------------------
p_check_buffer <- ggplot() +
  geom_sf(
    data = brazil_full_plot,
    fill = NA,
    color = "black",
    linewidth = 0.35
  ) +
  geom_sf(
    data = brazil_inner_plot,
    fill = NA,
    color = "red",
    linewidth = 0.35
  ) +
  coord_sf(
    xlim = xlim_br,
    ylim = ylim_br,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  theme_bw() +
  ggtitle("Check buffer: Brazil full boundary vs inner -25 km boundary")

ggsave(
  filename = file.path(dir_fig, "CHECK_Brazil_inner_buffer_25km.pdf"),
  plot = p_check_buffer,
  width = 12,
  height = 12,
  units = "cm"
)

# ==========================================================
# 9. RECORTAR/MASCARAR A BRASIL INTERIOR Y SUMAR 12 MESES
# ==========================================================
MB_Brazil_annual <- MB_Brazil_annual_raw |>
  crop(brazil_inner_vect) |>
  mask(brazil_inner_vect)

F60_Brazil_year <- stack60_year |>
  crop(brazil_inner_vect) |>
  mask(brazil_inner_vect)

F51_Brazil_year <- stack51_year |>
  crop(brazil_inner_vect) |>
  mask(brazil_inner_vect)

F60_Brazil_annual <- sum(F60_Brazil_year, na.rm = TRUE)
F51_Brazil_annual <- sum(F51_Brazil_year, na.rm = TRUE)

names(F60_Brazil_annual) <- paste0("MRBA60_", year_target)
names(F51_Brazil_annual) <- paste0("FireCCI51_", year_target)
names(MB_Brazil_annual)  <- paste0("MapBiomas_", year_target)

cat("\n--- SUMA TOTAL BRASIL INTERIOR", year_target, "(km²) ---\n")
cat("MapBiomas :", global(MB_Brazil_annual,  "sum", na.rm = TRUE)$sum, "\n")
cat("MRBA60    :", global(F60_Brazil_annual, "sum", na.rm = TRUE)$sum, "\n")
cat("FireCCI51 :", global(F51_Brazil_annual, "sum", na.rm = TRUE)$sum, "\n")

cat("\n--- GEOMETRÍAS DESPUÉS DEL RECORTE ---\n")
cat("MRBA60:\n")
print(F60_Brazil_annual)
cat("FireCCI51:\n")
print(F51_Brazil_annual)
cat("MapBiomas:\n")
print(MB_Brazil_annual)

# ==========================================================
# 9.3 ASEGURAR MISMA GEOMETRÍA PARA DIFERENCIAS Y CORRELACIÓN
# ==========================================================
if (!terra::compareGeom(F60_Brazil_annual, F51_Brazil_annual, stopOnError = FALSE)) {
  stop("MRBA60 y FireCCI51 no tienen la misma geometría tras el recorte.")
}

if (!terra::compareGeom(F60_Brazil_annual, MB_Brazil_annual, stopOnError = FALSE)) {
  
  warning(
    "MapBiomas no tiene exactamente la misma geometría que MRBA60 tras el recorte. ",
    "Se ajustará a la grilla de MRBA60 con resample(method='near'). ",
    "Revisa si esto es aceptable para tu validación."
  )
  
  MB_Brazil_annual <- terra::resample(
    MB_Brazil_annual,
    F60_Brazil_annual,
    method = "near"
  )
  
  MB_Brazil_annual <- mask(MB_Brazil_annual, F60_Brazil_annual)
  names(MB_Brazil_annual) <- paste0("MapBiomas_", year_target)
}

ext_br <- terra::ext(F60_Brazil_annual)

xlim_br_plot <- c(ext_br[1], ext_br[2])
ylim_br_plot <- c(ext_br[3], ext_br[4])

cat("\nDimensiones finales:\n")
cat("FireCCI51:", paste(dim(F51_Brazil_annual), collapse = " x "), "\n")
cat("MRBA60   :", paste(dim(F60_Brazil_annual), collapse = " x "), "\n")
cat("MapBiomas:", paste(dim(MB_Brazil_annual), collapse = " x "), "\n")

# ==========================================================
# 10. CORRELACIÓN SPEARMAN PIXEL A PIXEL
# ==========================================================
v51 <- as.vector(F51_Brazil_annual)
v60 <- as.vector(F60_Brazil_annual)
vMB <- as.vector(MB_Brazil_annual)

mask_cor <- !is.na(v51) & !is.na(v60) & !is.na(vMB)

# Opcional:
# mask_cor <- mask_cor & (v51 > 0 | v60 > 0 | vMB > 0)

v51 <- v51[mask_cor]
v60 <- v60[mask_cor]
vMB <- vMB[mask_cor]

cat("\nNúmero de píxeles usados en Spearman:", length(vMB), "\n")

if (length(vMB) < 3) {
  stop("Hay menos de 3 píxeles válidos para calcular Spearman.")
}

c51 <- cor.test(v51, vMB, method = "spearman", exact = FALSE)
c60 <- cor.test(v60, vMB, method = "spearman", exact = FALSE)

pvals <- c(c51$p.value, c60$p.value)
pvals_fdr <- p.adjust(pvals, method = "fdr")

spearman_cor_51 <- unname(c51$estimate)
spearman_cor_60 <- unname(c60$estimate)
p_value_FDR_51  <- pvals_fdr[1]
p_value_FDR_60  <- pvals_fdr[2]

cat("\n--- SPEARMAN VS MAPBIOMAS", year_target, "---\n")
cat("FireCCI51:", spearman_cor_51, " | p-FDR:", p_value_FDR_51, "\n")
cat("MRBA60   :", spearman_cor_60, " | p-FDR:", p_value_FDR_60, "\n")

# ==========================================================
# 11. PLOTS CON GGPLOT2: 4 MAPAS INDEPENDIENTES
# ==========================================================

# ----------------------------------------------------------
# 11.1 Copias para BA y diferencia
# ----------------------------------------------------------
MB_BA  <- MB_Brazil_annual
F60_BA <- F60_Brazil_annual
F51_BA <- F51_Brazil_annual

MB_BA[MB_BA == 0]   <- NA
F60_BA[F60_BA == 0] <- NA
F51_BA[F51_BA == 0] <- NA

MB_diff  <- MB_Brazil_annual
F60_diff <- F60_Brazil_annual

Diff_60_MB <- F60_diff - MB_diff
Diff_60_MB[Diff_60_MB == 0] <- NA

cat("\nRango diferencia MRBA60 - MapBiomas:\n")
print(terra::global(Diff_60_MB, "range", na.rm = TRUE))

# ----------------------------------------------------------
# 11.2 SpatRaster -> data.frame
# ----------------------------------------------------------
df_51 <- as.data.frame(F51_BA, xy = TRUE, na.rm = TRUE)
colnames(df_51)[3] <- "BA"

df_60 <- as.data.frame(F60_BA, xy = TRUE, na.rm = TRUE)
colnames(df_60)[3] <- "BA"

df_mb <- as.data.frame(MB_BA, xy = TRUE, na.rm = TRUE)
colnames(df_mb)[3] <- "BA"

df_diff <- as.data.frame(Diff_60_MB, xy = TRUE, na.rm = TRUE)
colnames(df_diff)[3] <- "diff_BA"

# ----------------------------------------------------------
# 11.3 Paleta y clases para BA
# ----------------------------------------------------------
cols_annual <- c(
  "#346099", "#5984BD", "#88B0E3", "#A9E5E8",
  "#FFFF8C", "#FDAE61", "#D73027", "#A50026", "#870000"
)

labels_annual <- c(
  "< 5", "5–10", "10–25", "25–50",
  "50–100", "100–250", "250–500", "500–750", "> 750"
)

breaks_annual <- c(0, 5, 10, 25, 50, 100, 250, 500, 750, Inf)

df_51 <- df_51 %>%
  mutate(
    BA_class = cut(
      BA,
      breaks = breaks_annual,
      labels = labels_annual,
      include.lowest = TRUE,
      right = FALSE
    )
  )

df_60 <- df_60 %>%
  mutate(
    BA_class = cut(
      BA,
      breaks = breaks_annual,
      labels = labels_annual,
      include.lowest = TRUE,
      right = FALSE
    )
  )

df_mb <- df_mb %>%
  mutate(
    BA_class = cut(
      BA,
      breaks = breaks_annual,
      labels = labels_annual,
      include.lowest = TRUE,
      right = FALSE
    )
  )

# ----------------------------------------------------------
# 11.4 Paleta continua para diferencia
# ----------------------------------------------------------
limit_diff <- 100

pal_diff <- colorRampPalette(c(
  "#08306B",
  "#2171B5",
  "#6BAED6",
  "#A8DADC",
  "#FFFAC0",
  "#FDAEFB",
  "#DD65B6",
  "#C51B7D",
  "#7A0177"
))(11)

vals <- seq(-100, 100, length.out = 11)

df_scale <- data.frame(
  x   = vals,
  col = pal_diff
)

pal_sv <- ggplot(df_scale) +
  geom_tile(aes(x = x, y = 1, fill = col)) +
  scale_fill_identity() +
  scale_x_continuous(breaks = c(-100, -50, 0, 50, 100)) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    text = element_text(size = 8)
  ) +
  xlab(expression(Delta~"Annual BA 2005 (km"^2*")"))

ggsave(
  filename = file.path(dir_fig, "palette_diff_Brazil_2005.pdf"),
  plot = pal_sv,
  width = 12,
  height = 2,
  units = "cm"
)

# ----------------------------------------------------------
# 11.5 Tema base común sin leyendas
# ----------------------------------------------------------
base_theme <- theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = 10),
    
    panel.grid.major = element_line(
      color = "grey85",
      linetype = "dotted",
      linewidth = 0.3
    ),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    
    axis.text.x = element_text(size = 9, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title  = element_text(size = 9, color = "black"),
    
    legend.position = "none",
    
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = ggplot2::margin(b = 2)
    ),
    
    # Márgenes reducidos para juntar más los paneles
    plot.margin = ggplot2::margin(2, 1, 2, 1)
  )

theme_no_y_labels <- theme(
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# ==========================================================
# 12. MAPAS
# ==========================================================

# ----------------------------------------------------------
# e) FireCCI51
# ----------------------------------------------------------
p_a <- ggplot() +
  geom_sf(
    data  = brazil_full_plot,
    fill  = NA,
    color = "black",
    linewidth = 0.25
  ) +
  geom_sf(
    data  = brazil_inner_plot,
    fill  = "grey20",
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_tile(
    data = df_51,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2005 (km²)",
    values = cols_annual,
    drop = FALSE,
    na.translate = FALSE
  ) +
  scale_x_continuous(
    name = "",
    breaks = lon_ticks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    breaks = lat_ticks,
    expand = c(0, 0)
  ) +
  coord_sf(
    xlim = xlim_br_plot,
    ylim = ylim_br_plot,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  ggtitle("e) FireCCI51")

# ----------------------------------------------------------
# f) MRBA60
# ----------------------------------------------------------
p_b <- ggplot() +
  geom_sf(
    data  = brazil_full_plot,
    fill  = NA,
    color = "black",
    linewidth = 0.25
  ) +
  geom_sf(
    data  = brazil_inner_plot,
    fill  = "grey20",
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_tile(
    data = df_60,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2005 (km²)",
    values = cols_annual,
    drop = FALSE,
    na.translate = FALSE
  ) +
  scale_x_continuous(
    name = "",
    breaks = lon_ticks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    breaks = lat_ticks,
    expand = c(0, 0)
  ) +
  coord_sf(
    xlim = xlim_br_plot,
    ylim = ylim_br_plot,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("f) MRBA60")

# ----------------------------------------------------------
# g) MapBiomas
# ----------------------------------------------------------
p_c <- ggplot() +
  geom_sf(
    data  = brazil_full_plot,
    fill  = NA,
    color = "black",
    linewidth = 0.25
  ) +
  geom_sf(
    data  = brazil_inner_plot,
    fill  = "grey20",
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_tile(
    data = df_mb,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2005 (km²)",
    values = cols_annual,
    drop = FALSE,
    na.translate = FALSE
  ) +
  scale_x_continuous(
    name = "",
    breaks = lon_ticks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    breaks = lat_ticks,
    expand = c(0, 0)
  ) +
  coord_sf(
    xlim = xlim_br_plot,
    ylim = ylim_br_plot,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("g) MapBiomas")

# ----------------------------------------------------------
# h) Diferencia MRBA60 - MapBiomas
# ----------------------------------------------------------
p_d <- ggplot() +
  geom_sf(
    data  = brazil_full_plot,
    fill  = NA,
    color = "black",
    linewidth = 0.25
  ) +
  geom_sf(
    data  = brazil_inner_plot,
    fill  = "grey20",
    color = "grey40",
    linewidth = 0.25
  ) +
  geom_tile(
    data = df_diff,
    aes(x = x, y = y, fill = diff_BA)
  ) +
  scale_fill_gradientn(
    colours = pal_diff,
    limits  = c(-limit_diff, limit_diff),
    oob     = scales::squish,
    name    = expression(Delta~"Annual BA 2005 (km"^2*")"~MRBA60 - MapBiomas)
  ) +
  scale_x_continuous(
    name = "",
    breaks = lon_ticks,
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name = "",
    breaks = lat_ticks,
    expand = c(0, 0)
  ) +
  coord_sf(
    xlim = xlim_br_plot,
    ylim = ylim_br_plot,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("h) MRBA60 - MapBiomas")

# ==========================================================
# 13. FIGURA 1×4 SIN LEYENDAS, MÁS COMPACTA
# ==========================================================
p_4panels <- p_a + p_b + p_c + p_d +
  plot_layout(
    ncol = 4,
    nrow = 1
  ) &
  theme(
    plot.margin = ggplot2::margin(2, 1, 2, 1)
  )

print(p_4panels)

ggsave(
  filename = file.path(
    dir_fig,
    "Figure5_Brazil_2005_4panels_1x4_noLegend.pdf"
  ),
  plot   = p_4panels,
  width  = 31,
  height = 8.5,
  units  = "cm"
)

ggsave(
  filename = file.path(
    dir_fig,
    "Figure5_Brazil_2005_4panels_1x4_noLegend.jpeg"
  ),
  plot   = p_4panels,
  width  = 31,
  height = 8.5,
  units  = "cm",
  dpi    = 600
)

# ==========================================================
# 14. GUARDAR CADA MAPA EN PDF Y JPEG INDEPENDIENTE
# ==========================================================
ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_e_F51.pdf"),
  plot     = p_a,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_e_F51.jpeg"),
  plot     = p_a,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_f_F60.pdf"),
  plot     = p_b,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_f_F60.jpeg"),
  plot     = p_b,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_g_MapBiomas.pdf"),
  plot     = p_c,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_g_MapBiomas.jpeg"),
  plot     = p_c,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_h_Diff60_MapBiomas.pdf"),
  plot     = p_d,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Brazil_2005_h_Diff60_MapBiomas.jpeg"),
  plot     = p_d,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

# ==========================================================
# 15. GUARDAR OBJETOS DE DIAGNÓSTICO
# ==========================================================
diag_buffer <- data.frame(
  year_target = year_target,
  buffer_dist_m = buffer_dist_m,
  crs_brazil_m = crs_brazil_m,
  area_brazil_full_km2 = area_brazil_full_km2,
  area_brazil_inner_km2 = area_brazil_inner_km2,
  area_removed_km2 = area_removed_km2,
  pct_area_removed = 100 * area_removed_km2 / area_brazil_full_km2,
  sum_MB_Brazil_inner_km2 = global(MB_Brazil_annual,  "sum", na.rm = TRUE)$sum,
  sum_MRBA60_Brazil_inner_km2 = global(F60_Brazil_annual, "sum", na.rm = TRUE)$sum,
  sum_FireCCI51_Brazil_inner_km2 = global(F51_Brazil_annual, "sum", na.rm = TRUE)$sum,
  spearman_FireCCI51_vs_MapBiomas = spearman_cor_51,
  spearman_MRBA60_vs_MapBiomas = spearman_cor_60,
  pFDR_FireCCI51_vs_MapBiomas = p_value_FDR_51,
  pFDR_MRBA60_vs_MapBiomas = p_value_FDR_60,
  n_pixels_spearman = length(vMB)
)

write.csv(
  diag_buffer,
  file = file.path(dir_fig, "CHECK_Brazil_2005_buffer25km_diagnostics.csv"),
  row.names = FALSE
)

save(
  brazil_sf_full,
  brazil_inner_sf,
  MB_Brazil_annual,
  F60_Brazil_annual,
  F51_Brazil_annual,
  Diff_60_MB,
  diag_buffer,
  file = file.path(dir_fig, "CHECK_Brazil_2005_buffer25km_objects.RData")
)

cat("\nProceso finalizado correctamente.\n")
cat("Figuras y diagnósticos guardados en:\n")
cat(dir_fig, "\n")