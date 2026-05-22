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

# ==========================================================
# 1. CARGAR S2 2016 (burned_area en m²) Y PASAR A km²
# ==========================================================
path_s2 <- "/mnt/disco6tb/Validacion/AF_S2_2016"

files_s2 <- sort(list.files(path_s2, pattern = "\\.nc$", full.names = TRUE))
files_s2

if (length(files_s2) == 0) {
  stop("No se encontraron archivos NetCDF de S2 en: ", path_s2)
}

# Leemos sólo la variable burned_area (m²) de cada mes
s2_list_m2 <- lapply(files_s2, function(f) {
  rast(f, sub = "burned_area")
})

# Stack mensual S2 2016 (12 capas, m²)
s2_stack_m2 <- rast(s2_list_m2)

# Convertir TODO S2 a km²
s2_stack_km2 <- s2_stack_m2 / 1e6
names(s2_stack_km2) <- sprintf("S2_2016_%02d", 1:nlyr(s2_stack_km2))

# Plantilla geométrica definida por S2
template <- s2_stack_km2[[1]]

cat("\nTemplate S2:\n")
print(template)

# Suma global S2 2016
s2_2016_global <- global(
  sum(s2_stack_km2, na.rm = TRUE),
  "sum",
  na.rm = TRUE
)$sum

cat("Suma global S2 2016 (km²):", s2_2016_global, "\n")

# ==========================================================
# 2. CARGAR LAT/LON Y ARRAYS MRBA60 / FireCCI51 (km²)
# ==========================================================
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")

cat("\nRango lat:", paste(range(lat), collapse = " - "), "\n")
cat("Rango lon:", paste(range(lon), collapse = " - "), "\n")

Modelo <- "B1-MRBA60-2003-2024"

output_base <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
output_dir_RData <- file.path(output_base, "RData")

# ----------------------------------------------------------
# MRBA60 armonizado
# ----------------------------------------------------------
ruta_RData <- file.path(
  output_dir_RData,
  "MRBA60_BA_m2_monthly_2003_2024.RData"
)

load(ruta_RData)

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

cat("length(lon):", length(lon), " | ncol Fire arrays:", dim(BA_MRBA60)[1], "\n")
cat("length(lat):", length(lat), " | nrow Fire arrays:", dim(BA_MRBA60)[2], "\n")

stopifnot(dim(BA_MRBA60)[1] == length(lon))
stopifnot(dim(BA_MRBA60)[2] == length(lat))
stopifnot(all(dim(BA_MRBA60) == dim(BA_FireCCI51)))

# ==========================================================
# 3. ÍNDICES DE LOS 12 MESES DE 2016 EN LOS ARRAYS
#    Serie mensual desde 2003-01
# ==========================================================
start_year  <- 2003
year_target <- 2016

start_idx_2016 <- (year_target - start_year) * 12 + 1
end_idx_2016   <- start_idx_2016 + 11

idx_2016 <- start_idx_2016:end_idx_2016

cat("\nÍndices 2016:", paste(idx_2016, collapse = ","), "\n")

arr60_2016 <- BA_MRBA60[, , idx_2016, drop = FALSE]
arr51_2016 <- BA_FireCCI51[, , idx_2016, drop = FALSE]

cat("Dim arr60_2016:", paste(dim(arr60_2016), collapse = " x "), "\n")
cat("Dim arr51_2016:", paste(dim(arr51_2016), collapse = " x "), "\n")

cat("Rango arr60_2016:", paste(range(arr60_2016, na.rm = TRUE), collapse = " - "), "\n")
cat("Rango arr51_2016:", paste(range(arr51_2016, na.rm = TRUE), collapse = " - "), "\n")

sum_arr60_2016 <- sum(arr60_2016, na.rm = TRUE)
sum_arr51_2016 <- sum(arr51_2016, na.rm = TRUE)

cat("Suma global array MRBA60 2016 (km²):", sum_arr60_2016, "\n")
cat("Suma global array FireCCI51 2016 (km²):", sum_arr51_2016, "\n")

# ==========================================================
# 4. RASTERIZAR MRBA60 2016 CON LA GEOMETRÍA DE S2
# ==========================================================
ncols    <- ncol(template)
nrows    <- nrow(template)
ext_glob <- ext(template)
crs_glob <- crs(template)

cat("\ntemplate ncols/nrows:", ncols, nrows, "\n")
cat("length(lon)/length(lat):", length(lon), length(lat), "\n")

stopifnot(ncols == length(lon))
stopifnot(nrows == length(lat))

stack60_2016 <- rast(
  ncols   = ncols,
  nrows   = nrows,
  nlyrs   = 12,
  extent  = ext_glob,
  crs     = crs_glob
)

for (m in 1:12) {
  cat("Rellenando MRBA60 mes", m, "\n")
  
  slice_2d <- arr60_2016[, , m]        # [lon, lat]
  mat_latlon <- t(slice_2d)            # [lat, lon]
  mat_latlon_flip <- mat_latlon[nrow(mat_latlon):1, ]
  
  stack60_2016[[m]][] <- as.vector(t(mat_latlon_flip))
}

names(stack60_2016) <- sprintf("MRBA60_2016_%02d", 1:12)

sum_rast60_2016 <- global(
  sum(stack60_2016, na.rm = TRUE),
  "sum",
  na.rm = TRUE
)$sum

cat("Suma global raster MRBA60 2016 (km²):", sum_rast60_2016, "\n")

# ==========================================================
# 5. RASTERIZAR FireCCI51 2016 CON LA MISMA GEOMETRÍA
# ==========================================================
stack51_2016 <- rast(
  ncols   = ncols,
  nrows   = nrows,
  nlyrs   = 12,
  extent  = ext_glob,
  crs     = crs_glob
)

for (m in 1:12) {
  cat("Rellenando FireCCI51 mes", m, "\n")
  
  slice_2d <- arr51_2016[, , m]        # [lon, lat]
  mat_latlon <- t(slice_2d)            # [lat, lon]
  mat_latlon_flip <- mat_latlon[nrow(mat_latlon):1, ]
  
  stack51_2016[[m]][] <- as.vector(t(mat_latlon_flip))
}

names(stack51_2016) <- sprintf("FireCCI51_2016_%02d", 1:12)

sum_rast51_2016 <- global(
  sum(stack51_2016, na.rm = TRUE),
  "sum",
  na.rm = TRUE
)$sum

cat("Suma global raster FireCCI51 2016 (km²):", sum_rast51_2016, "\n")

# ==========================================================
# 6. RECORTAR ÁFRICA A LAT [-35, 20] Y SUMAR LOS 12 MESES
# ==========================================================

africa_sf_full <- ne_countries(
  continent = "Africa",
  returnclass = "sf"
)

africa_sf_full <- africa_sf_full |>
  st_transform(4326) |>
  st_make_valid()

africa_vect_full <- vect(africa_sf_full)

lat_min <- -35
lat_max <-  20

ext_africa_full <- terra::ext(africa_vect_full)

ext_recorte <- terra::ext(
  ext_africa_full[1],
  ext_africa_full[2],
  lat_min,
  lat_max
)

africa_vect <- terra::crop(africa_vect_full, ext_recorte)
africa_sf <- sf::st_as_sf(africa_vect)

S2_Africa_2016  <- mask(crop(s2_stack_km2, africa_vect), africa_vect)
F60_Africa_2016 <- mask(crop(stack60_2016, africa_vect), africa_vect)
F51_Africa_2016 <- mask(crop(stack51_2016, africa_vect), africa_vect)

S2_Africa_annual  <- sum(S2_Africa_2016, na.rm = TRUE)
F60_Africa_annual <- sum(F60_Africa_2016, na.rm = TRUE)
F51_Africa_annual <- sum(F51_Africa_2016, na.rm = TRUE)

cat("\n--- SUMA TOTAL ÁFRICA 2016 (km²) EN LAT [-35, 20] ---\n")
cat("S2   :", global(S2_Africa_annual, "sum", na.rm = TRUE)$sum, "\n")
cat("F60  :", global(F60_Africa_annual, "sum", na.rm = TRUE)$sum, "\n")
cat("F51  :", global(F51_Africa_annual, "sum", na.rm = TRUE)$sum, "\n")

ext_af  <- terra::ext(S2_Africa_annual)
xlim_af <- c(ext_af[1], ext_af[2])
ylim_af <- c(ext_af[3], ext_af[4])

# Longitud de 10 en 10
lon_ticks <- seq(-20, 60, 10)

# Latitud de 5 en 5
lat_ticks <- seq(-35, 20, 5)

cat("\nDimensiones finales:\n")
cat("FireCCI51:", paste(dim(F51_Africa_annual), collapse = " x "), "\n")
cat("MRBA60  :", paste(dim(F60_Africa_annual), collapse = " x "), "\n")
cat("S2       :", paste(dim(S2_Africa_annual), collapse = " x "), "\n")

# ==========================================================
# 6.1 CORRELACIÓN SPEARMAN PIXEL A PIXEL
# ==========================================================
v51 <- as.vector(F51_Africa_annual)
v60 <- as.vector(F60_Africa_annual)
vS2 <- as.vector(S2_Africa_annual)

mask_cor <- !is.na(v51) & !is.na(v60) & !is.na(vS2)

v51 <- v51[mask_cor]
v60 <- v60[mask_cor]
vS2 <- vS2[mask_cor]

cat("\nNúmero de píxeles usados en Spearman:", length(vS2), "\n")

if (length(vS2) < 3) {
  stop("Hay menos de 3 píxeles válidos para calcular Spearman.")
}

c51 <- cor.test(v51, vS2, method = "spearman", exact = FALSE)
c60 <- cor.test(v60, vS2, method = "spearman", exact = FALSE)

pvals <- c(c51$p.value, c60$p.value)
pvals_fdr <- p.adjust(pvals, method = "fdr")

spearman_cor_51 <- unname(c51$estimate)
spearman_cor_60 <- unname(c60$estimate)
p_value_FDR_51  <- pvals_fdr[1]
p_value_FDR_60  <- pvals_fdr[2]

cat("\n--- SPEARMAN VS FireCCISFD11 2016 ---\n")
cat("FireCCI51:", spearman_cor_51, " | p-FDR:", p_value_FDR_51, "\n")
cat("MRBA60  :", spearman_cor_60, " | p-FDR:", p_value_FDR_60, "\n")

# ==========================================================
# 7. PLOTS CON GGPLOT2: 4 MAPAS INDEPENDIENTES a)–d)
# ==========================================================

# ----------------------------------------------------------
# 7.1 COPIAS PARA BA Y DIFERENCIA
# ----------------------------------------------------------
S2_BA  <- S2_Africa_annual
F60_BA <- F60_Africa_annual
F51_BA <- F51_Africa_annual

S2_BA[S2_BA == 0]   <- NA
F60_BA[F60_BA == 0] <- NA
F51_BA[F51_BA == 0] <- NA

S2_diff  <- S2_Africa_annual
F60_diff <- F60_Africa_annual

Diff_60_S2 <- F60_diff - S2_diff
Diff_60_S2[Diff_60_S2 == 0] <- NA

cat("\nRango diferencia MRBA60 - FireCCISFD11:\n")
print(terra::global(Diff_60_S2, "range", na.rm = TRUE))

# ----------------------------------------------------------
# 7.2 SpatRaster -> data.frame
# ----------------------------------------------------------
df_51 <- as.data.frame(F51_BA, xy = TRUE, na.rm = TRUE)
colnames(df_51)[3] <- "BA"

df_60 <- as.data.frame(F60_BA, xy = TRUE, na.rm = TRUE)
colnames(df_60)[3] <- "BA"

df_s2 <- as.data.frame(S2_BA, xy = TRUE, na.rm = TRUE)
colnames(df_s2)[3] <- "BA"

df_diff <- as.data.frame(Diff_60_S2, xy = TRUE, na.rm = TRUE)
colnames(df_diff)[3] <- "diff_BA"

# ----------------------------------------------------------
# 7.3 PALETA Y CLASES PARA BA
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

df_s2 <- df_s2 %>%
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
# 7.4 PALETA CONTINUA PARA DIFERENCIA
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

dir_fig <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure5_R1"
dir.create(dir_fig, recursive = TRUE, showWarnings = FALSE)

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
  xlab(expression(Delta~"Annual BA 2016 (km"^2*")"))

ggsave(
  filename = file.path(dir_fig, "palette_diff_Africa_2016.pdf"),
  plot = pal_sv,
  width = 12,
  height = 2,
  units = "cm"
)

# ----------------------------------------------------------
# 7.5 TEMA BASE COMÚN
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
    
    plot.margin = ggplot2::margin(2, 1, 2, 1)
  )

theme_no_y_labels <- theme(
  axis.text.y  = element_blank(),
  axis.ticks.y = element_blank()
)

# ==========================================================
# 8. MAPAS a)–d)
# ==========================================================

# ----------------------------------------------------------
# a) FireCCI51
# ----------------------------------------------------------
p_a <- ggplot() +
  geom_sf(
    data  = africa_sf,
    fill  = "grey20",
    color = "grey0",
    linewidth = 0.3
  ) +
  geom_tile(
    data = df_51,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2016 (km²)",
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
    xlim = xlim_af,
    ylim = ylim_af,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  ggtitle("a) FireCCI51")

# ----------------------------------------------------------
# b) MRBA60
# ----------------------------------------------------------
p_b <- ggplot() +
  geom_sf(
    data  = africa_sf,
    fill  = "grey20",
    color = "grey0",
    linewidth = 0.3
  ) +
  geom_tile(
    data = df_60,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2016 (km²)",
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
    xlim = xlim_af,
    ylim = ylim_af,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("b) MRBA60")

# ----------------------------------------------------------
# c) FireCCISFD11
# ----------------------------------------------------------
p_c <- ggplot() +
  geom_sf(
    data  = africa_sf,
    fill  = "grey20",
    color = "grey0",
    linewidth = 0.3
  ) +
  geom_tile(
    data = df_s2,
    aes(x = x, y = y, fill = BA_class)
  ) +
  scale_fill_manual(
    name = "Annual BA 2016 (km²)",
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
    xlim = xlim_af,
    ylim = ylim_af,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("c) FireCCISFD11")

# ----------------------------------------------------------
# d) MRBA60 - FireCCISFD11
# ----------------------------------------------------------
p_d <- ggplot() +
  geom_sf(
    data  = africa_sf,
    fill  = "grey20",
    color = "grey0",
    linewidth = 0.3
  ) +
  geom_tile(
    data = df_diff,
    aes(x = x, y = y, fill = diff_BA)
  ) +
  scale_fill_gradientn(
    colours = pal_diff,
    limits  = c(-limit_diff, limit_diff),
    oob     = scales::squish,
    name    = expression(Delta~"Annual BA 2016 (km"^2*")"~MRBA60 - FireCCISFD11)
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
    xlim = xlim_af,
    ylim = ylim_af,
    expand = FALSE,
    crs = sf::st_crs(4326)
  ) +
  base_theme +
  theme_no_y_labels +
  ggtitle("d) MRBA60 - FireCCISFD11")

# ==========================================================
# 9. FIGURA 1×4 SIN LEYENDAS, MÁS COMPACTA
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
    "Figure5_Africa_4panels_1x4_noLegend.pdf"
  ),
  plot   = p_4panels,
  width  = 31,
  height = 8.5,
  units  = "cm"
)

ggsave(
  filename = file.path(
    dir_fig,
    "Figure5_Africa_4panels_1x4_noLegend.jpeg"
  ),
  plot   = p_4panels,
  width  = 31,
  height = 8.5,
  units  = "cm",
  dpi    = 600
)

# ==========================================================
# 10. GUARDAR CADA MAPA EN PDF Y JPEG INDEPENDIENTE
# ==========================================================
ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_a_F51.pdf"),
  plot     = p_a,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_a_F51.jpeg"),
  plot     = p_a,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_b_F60.pdf"),
  plot     = p_b,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_b_F60.jpeg"),
  plot     = p_b,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_c_S2.pdf"),
  plot     = p_c,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_c_S2.jpeg"),
  plot     = p_c,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_d_Diff60_S2.pdf"),
  plot     = p_d,
  width    = 8.5,
  height   = 8.5,
  units    = "cm"
)

ggsave(
  filename = file.path(dir_fig, "Figure5_Africa_d_Diff60_S2.jpeg"),
  plot     = p_d,
  width    = 8.5,
  height   = 8.5,
  units    = "cm",
  dpi      = 600
)

# ==========================================================
# 11. GUARDAR DIAGNÓSTICOS
# ==========================================================
diag_africa <- data.frame(
  year_target = year_target,
  lat_min = lat_min,
  lat_max = lat_max,
  sum_S2_Africa_km2 = global(S2_Africa_annual, "sum", na.rm = TRUE)$sum,
  sum_MRBA60_Africa_km2 = global(F60_Africa_annual, "sum", na.rm = TRUE)$sum,
  sum_FireCCI51_Africa_km2 = global(F51_Africa_annual, "sum", na.rm = TRUE)$sum,
  spearman_FireCCI51_vs_FireCCISFD11 = spearman_cor_51,
  spearman_MRBA60_vs_FireCCISFD11 = spearman_cor_60,
  pFDR_FireCCI51_vs_FireCCISFD11 = p_value_FDR_51,
  pFDR_MRBA60_vs_FireCCISFD11 = p_value_FDR_60,
  n_pixels_spearman = length(vS2)
)

write.csv(
  diag_africa,
  file = file.path(dir_fig, "CHECK_Africa_2016_diagnostics.csv"),
  row.names = FALSE
)

save(
  S2_Africa_annual,
  F60_Africa_annual,
  F51_Africa_annual,
  Diff_60_S2,
  diag_africa,
  file = file.path(dir_fig, "CHECK_Africa_2016_objects.RData")
)

cat("\nProceso finalizado correctamente.\n")
cat("Figuras y diagnósticos guardados en:\n")
cat(dir_fig, "\n")