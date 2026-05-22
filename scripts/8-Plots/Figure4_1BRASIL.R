
# ==========================================================
# # LIMPIEZA INICIAL
# # ==========================================================
rm(list = ls())
graphics.off()
gc()

# ==========================================================
# LIBRERÍAS
# ==========================================================
libs <- c("terra", "ncdf4", "dplyr", "ggplot2", "tidyr",
          "purrr", "tibble", "scales", "raster", "sf",
          "rnaturalearth", "rnaturalearthdata")

for (p in libs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}

sf::sf_use_s2(FALSE)

# ==========================================================
# CONFIGURACIÓN DE RUTAS Y MODELO
# ==========================================================
Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
output_base <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv   <- file.path(output_base, "csv")
output_dir_plot  <- file.path(output_base, "plot")
output_dir_RData <- file.path(output_base, "RData")

dir.create(output_base,      showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_csv,   showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_plot,  showWarnings = FALSE, recursive = TRUE)
dir.create(output_dir_RData, showWarnings = FALSE, recursive = TRUE)

dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure4_R1/"
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

col_map <- c(
  "MRBA60" = "blue",
  "FireCCI51" = "#E69F00",
  "MCD64A1"   = "brown4",
  "GFED5"     = "gray30"
)

# ==========================================================
# HELPERS
# ==========================================================
make_years_mensual <- function(start_year, ntime) {
  stopifnot(ntime %% 12 == 0)
  rep(start_year:(start_year + ntime/12 - 1), each = 12)
}

sumar_mensual_a_anual <- function(arr_mensual, years_mensual, years_unicos) {
  stopifnot(length(dim(arr_mensual)) == 3)
  nx <- dim(arr_mensual)[1]; ny <- dim(arr_mensual)[2]
  out <- array(NA_real_, dim = c(nx, ny, length(years_unicos)))
  for (j in seq_along(years_unicos)) {
    idx <- which(years_mensual == years_unicos[j])
    out[, , j] <- apply(arr_mensual[, , idx, drop = FALSE], c(1, 2), sum, na.rm = TRUE)
  }
  out
}
# ==========================================================
# Cargar armonizado (MRBA60) mensual 0.25º
# ==========================================================
ruta_RData <- file.path( output_dir_RData, "MRBA60_BA_m2_monthly_2003_2024.RData" )
load(ruta_RData)
BA_final <- BA_MRBA60/ 1e6

# BA_final <- BA_FireCCI60 / 1e6  # m2 -> km²
BA_final[is.na(BA_final)] <- 0

years_mensual_fc60 <- make_years_mensual(2003, dim(BA_final)[3])
years_unicos_fc60  <- unique(years_mensual_fc60)

# ==========================================================
# Cargar lon/lat (0.25º)
# ==========================================================
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")
lon_glob <- lon
lat_glob <- lat

# ==========================================================
# Cargar productos BA externos (mensuales 0.25º)
# ==========================================================
# FireCCI51 (2001–2022) -> submuestreamos 2003–2022
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))
BA_FireCCI51 <- f51 / 1e6
BA_FireCCI51[is.na(BA_FireCCI51)] <- 0

years_mensual_fc51 <- make_years_mensual(2003, dim(BA_FireCCI51)[3])
years_unicos_fc51  <- unique(years_mensual_fc51)


nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/MCD64CMQ_Monthly_2000-2024.nc"
nc_mcd <- nc_open(nc_path_mcd)
BA_MCD64 <- ncvar_get(nc_mcd, "band_data") / 1e6  # (según tu comentario previo)
nc_close(nc_mcd)
BA_MCD64 <- BA_MCD64[ , , 27:290, drop = FALSE]  # 2003–2024 (264)
BA_MCD64[is.na(BA_MCD64)] <- 0
dim(BA_MCD64)

years_mensual_mcd64 <- make_years_mensual(2003, dim(BA_MCD64)[3])
years_unicos_mcd64  <- unique(years_mensual_mcd64)
# GFED5 (2003–2024 en tu netcdf)
nc_gfed <- nc_open("/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc")
BA_GFED5 <- ncvar_get(nc_gfed, "burned_area") / 1e6  # ajusta si tu variable ya está en km²
nc_close(nc_gfed)
BA_GFED5[is.na(BA_GFED5)] <- 0

years_mensual_gfed5 <- make_years_mensual(2003, dim(BA_GFED5)[3])
years_unicos_gfed5  <- unique(years_mensual_gfed5)

# # ==========================================================
# # CARGAR ARMONIZADO (FireCCI60) MENSUAL 0.25º
# # ==========================================================
# load("/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-3/RData/FireCCI60_BA_m2_monthly_2003_2024.RData")
# BA_final <- BA_FireCCI60 / 1e6  # m2 -> km²
# BA_final[is.na(BA_final)] <- 0
# 
# years_mensual_fc60 <- make_years_mensual(2003, dim(BA_final)[3])
# years_unicos_fc60  <- unique(years_mensual_fc60)
# 
# # ==========================================================
# # CARGAR lon/lat (0.25º)
# # ==========================================================
# load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/longitude.RData")
# load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/latitude.RData")
# lon_glob <- lon
# lat_glob <- lat
# 
# # ==========================================================
# # CARGAR PRODUCTOS BA EXTERNOS (MENSUALES 0.25º)
# # ==========================================================
# # FireCCI51 (2001–2022) -> submuestreamos 2003–2022
# load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/FireCCI51_2001_2022_0.25degree-download.RData")
# BA_FireCCI51 <- Fire51 / 1e6
# BA_FireCCI51 <- BA_FireCCI51[,,25:264]  # 2003–2022
# BA_FireCCI51[is.na(BA_FireCCI51)] <- 0
# years_mensual_fc51 <- make_years_mensual(2003, dim(BA_FireCCI51)[3])
# years_unicos_fc51  <- unique(years_mensual_fc51)
# 
# # MCD64A1 (2001–2022) -> submuestreamos 2003–2022
# nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/nc_output/MCD64CMQ_BA_2001_2022.nc"
# nc_mcd <- nc_open(nc_path_mcd)
# BA_MCD64 <- ncvar_get(nc_mcd, "Band1") / 1e4
# BA_MCD64 <- BA_MCD64[,,25:264]  # 2003–2022
# nc_close(nc_mcd)
# BA_MCD64[is.na(BA_MCD64)] <- 0
# years_mensual_mcd64 <- make_years_mensual(2003, dim(BA_MCD64)[3])
# years_unicos_mcd64  <- unique(years_mensual_mcd64)
# 
# # GFED5 (2003–2024 según tu nc)
# nc_gfed <- nc_open("/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc")
# BA_GFED5 <- ncvar_get(nc_gfed, "burned_area") / 1e6  # ajusta si ya está en km²
# nc_close(nc_gfed)
# BA_GFED5[is.na(BA_GFED5)] <- 0
# years_mensual_gfed5 <- make_years_mensual(2003, dim(BA_GFED5)[3])
# years_unicos_gfed5  <- unique(years_mensual_gfed5)

# ==========================================================
# PRECOMPUTAR ANUAL PARA LOS PRODUCTOS GLOBALES (cada uno con su tiempo)
# ==========================================================
BA_fc60_anual  <- sumar_mensual_a_anual(BA_final,     years_mensual_fc60,  years_unicos_fc60)
BA_fc51_anual  <- sumar_mensual_a_anual(BA_FireCCI51, years_mensual_fc51,  years_unicos_fc51)
BA_mcd64_anual <- sumar_mensual_a_anual(BA_MCD64,     years_mensual_mcd64, years_unicos_mcd64)
BA_gfed5_anual <- sumar_mensual_a_anual(BA_GFED5,     years_mensual_gfed5, years_unicos_gfed5)

# ==========================================================
# MÁSCARA DE BRASIL INTERIOR EN REJILLA GLOBAL (0.25º)
# ==========================================================
brasil_sf <- rnaturalearth::ne_countries(country = "Brazil", scale = "medium", returnclass = "sf")
brasil_proj  <- st_transform(brasil_sf, 3857)
brasil_union <- st_union(brasil_proj)

buffer_dist_m <- 25000  # 25 km hacia dentro
brasil_inner_proj <- st_buffer(brasil_union, dist = -buffer_dist_m)
brasil_inner <- st_transform(brasil_inner_proj, 4326)

ix_lon_sa <- which(lon_glob >= -90 & lon_glob <= -20)
ix_lat_sa <- which(lat_glob >= -40 & lat_glob <= 15)

lon_sa <- lon_glob[ix_lon_sa]
lat_sa <- lat_glob[ix_lat_sa]

grid_sa <- expand.grid(lon = lon_sa, lat = lat_sa)
pts_sa  <- st_as_sf(grid_sa, coords = c("lon","lat"), crs = 4326)

inside_sa <- st_within(pts_sa, brasil_inner, sparse = FALSE)[,1]

mask_sa <- matrix(inside_sa,
                  nrow = length(lon_sa),
                  ncol = length(lat_sa),
                  byrow = FALSE)

mask_brasil <- matrix(FALSE, nrow = length(lon_glob), ncol = length(lat_glob))
mask_brasil[ix_lon_sa, ix_lat_sa] <- mask_sa

save(mask_brasil, file = file.path(output_dir_RData, "mask_brasil_025deg_inner.RData"))

# ==========================================================
# VALIDACIÓN MapBiomas (BRASIL, 0.25º, km²/año) + MÁSCARA INTERIOR
# ==========================================================
ruta_MapBrasil <- "/mnt/disco6tb/Validacion/MapBrasil"
patron_mb   <- "^burned_km2_025deg_(\\d{4})\\.tif$"
archivos_mb <- list.files(path = ruta_MapBrasil, pattern = patron_mb, full.names = TRUE)

if (length(archivos_mb) == 0) stop("No se encontraron archivos de MapBiomas en: ", ruta_MapBrasil)
archivos_mb <- archivos_mb[order(archivos_mb)]
years_mb <- as.numeric(sub(patron_mb, "\\1", basename(archivos_mb)))

stk_mb <- stack(archivos_mb)
e_mb   <- extent(stk_mb)
dim_x  <- ncol(stk_mb)
dim_y  <- nrow(stk_mb)
nt_mb  <- nlayers(stk_mb)

# Recorte global a bbox MapBiomas
ix_lon <- which(lon_glob >= e_mb@xmin & lon_glob <= e_mb@xmax)
ix_lat <- which(lat_glob >= e_mb@ymin & lat_glob <= e_mb@ymax)

lon_br <- lon_glob[ix_lon]
lat_br <- lat_glob[ix_lat]

mask_br_subset <- mask_brasil[ix_lon, ix_lat]

# Array MapBiomas [lon, lat, time] (ya viene en km²/año según patrón de archivo)
BA_MapBrasil <- array(NA_real_, dim = c(length(lon_br), length(lat_br), nt_mb))

for (i in seq_len(nt_mb)) {
  mat <- raster::as.matrix(stk_mb[[i]], na.rm = FALSE)
  
  if (!is.matrix(mat)) {
    mat <- matrix(mat, nrow = dim_y, ncol = dim_x, byrow = TRUE)
  } else if (nrow(mat) != dim_y || ncol(mat) != dim_x) {
    mat <- matrix(as.vector(mat), nrow = dim_y, ncol = dim_x, byrow = TRUE)
  }
  dimnames(mat) <- NULL
  
  mat_flip <- mat[dim_y:1, , drop = FALSE]
  slice    <- t(mat_flip)
  dimnames(slice) <- NULL
  
  if (nrow(slice) != length(lon_br) || ncol(slice) != length(lat_br)) {
    stop("Dimensiones MapBiomas no coinciden con lon_br/lat_br. Revisa rejilla.")
  }
  
  BA_MapBrasil[,,i] <- slice
}

# Aplicar máscara Brasil interior
for (i in seq_len(nt_mb)) {
  tmp <- BA_MapBrasil[,,i]
  tmp[!mask_br_subset] <- NA
  BA_MapBrasil[,,i] <- tmp
}

# Recorte + máscara de productos globales (anuales)
fc60_br  <- BA_fc60_anual[ix_lon, ix_lat, , drop = FALSE]
fc51_br  <- BA_fc51_anual[ix_lon, ix_lat, , drop = FALSE]
mcd64_br <- BA_mcd64_anual[ix_lon, ix_lat, , drop = FALSE]
gfed5_br <- BA_gfed5_anual[ix_lon, ix_lat, , drop = FALSE]

for (t in seq_len(dim(fc51_br)[3])) {
  fc60_br[,,t][!mask_br_subset]  <- NA
  fc51_br[,,t][!mask_br_subset]  <- NA
  mcd64_br[,,t][!mask_br_subset] <- NA
  gfed5_br[,,t][!mask_br_subset] <- NA
}

length(mask_br_subset)

# ==========================================================
# AÑOS COMUNES (FIJAR 2003–2022 Y VERIFICAR QUE DE VERDAD EXISTEN)
# ==========================================================
start_year_common <- 2003
end_year_common   <- 2024

anos_comunes <- Reduce(intersect, list(
  years_mb,
  years_unicos_fc60,
  years_unicos_fc51,
  years_unicos_mcd64,
  years_unicos_gfed5
))
anos_comunes <- sort(anos_comunes)
anos_comunes <- anos_comunes[anos_comunes >= start_year_common & anos_comunes <= end_year_common]

cat("MapBiomas years range: ", min(years_mb), "-", max(years_mb), "\n")
cat("Common years used range: ", min(anos_comunes), "-", max(anos_comunes), "\n")
cat("N common years: ", length(anos_comunes), "\n")
cat("First years: ", paste(head(anos_comunes, 5), collapse = ", "), "\n")
cat("Last years:  ", paste(tail(anos_comunes, 5), collapse = ", "), "\n")

stopifnot(length(anos_comunes) > 0)
stopifnot(min(anos_comunes) == start_year_common)
stopifnot(max(anos_comunes) == end_year_common)
stopifnot(length(anos_comunes) == (end_year_common - start_year_common + 1))

idx_map   <- match(anos_comunes, years_mb)
idx_fc60  <- match(anos_comunes, years_unicos_fc60)
idx_fc51  <- match(anos_comunes, years_unicos_fc51)
idx_mcd64 <- match(anos_comunes, years_unicos_mcd64)
idx_gfed5 <- match(anos_comunes, years_unicos_gfed5)

# ==========================================================
# SUMAS ANUALES (SOLO BRASIL INTERIOR)
# ==========================================================
nA <- length(anos_comunes)
sum_map   <- sum_fc60 <- sum_fc51 <- sum_mcd64 <- sum_gfed5 <- numeric(nA)

for (k in seq_len(nA)) {
  im  <- idx_map[k]
  i60 <- idx_fc60[k]
  i51 <- idx_fc51[k]
  imc <- idx_mcd64[k]
  igf <- idx_gfed5[k]
  
  sum_map[k]   <- sum(BA_MapBrasil[,,im], na.rm = TRUE)
  sum_fc60[k]  <- sum(fc60_br[,,i60],     na.rm = TRUE)
  sum_fc51[k]  <- sum(fc51_br[,,i51],     na.rm = TRUE)
  sum_mcd64[k] <- sum(mcd64_br[,,imc],    na.rm = TRUE)
  sum_gfed5[k] <- sum(gfed5_br[,,igf],    na.rm = TRUE)
}

# ==========================================================
# DATA FRAME
# ==========================================================
df_mb <- data.frame(
  year  = anos_comunes,
  local = sum_map,
  fc60  = sum_fc60,
  fc51  = sum_fc51,
  mcd64 = sum_mcd64,
  gfed5 = sum_gfed5
)

df_mb_long <- df_mb %>%
  pivot_longer(cols = c(fc60, fc51, mcd64, gfed5),
               names_to = "dataset", values_to = "value") %>%
  mutate(
    dataset = factor(dataset,
                     levels = c("fc60", "fc51", "mcd64", "gfed5"),
                     labels = c("MRBA60", "FireCCI51", "MCD64A1", "GFED5")),
    region = "Brasil (MapBiomas, interior)"
  )

# ==========================================================
# DIAGNÓSTICO: ¿POR QUÉ MAE Y BIAS PODRÍAN SER IGUALES?
# ==========================================================
diag_err <- df_mb_long %>%
  mutate(err = value - local) %>%
  group_by(dataset) %>%
  summarise(
    n_years   = n(),
    min_year  = min(year),
    max_year  = max(year),
    min_err   = min(err, na.rm = TRUE),
    max_err   = max(err, na.rm = TRUE),
    n_neg     = sum(err < 0, na.rm = TRUE),
    n_pos     = sum(err > 0, na.rm = TRUE),
    n_zero    = sum(err == 0, na.rm = TRUE),
    mae_raw   = mean(abs(err), na.rm = TRUE),
    bias_raw  = mean(err, na.rm = TRUE),
    .groups = "drop"
  )

print(diag_err)

write.csv(diag_err,
          file = file.path(dir_plot,
                           paste0("diagnostico_err_MapBiomas_Brasil_",
                                  start_year_common, "_", end_year_common, ".csv")),
          row.names = FALSE)

# ==========================================================
# MÉTRICAS (incluye BIAS % + Spearman con significancia 95%)
# ==========================================================
groups_mb <- df_mb_long %>% group_split(dataset, .keep = TRUE)

metrics_mb <- purrr::map_dfr(groups_mb, function(d){
  ds <- as.character(unique(d$dataset))
  if (nrow(d) == 0) {
    return(tibble(
      dataset      = ds,
      r2           = NA_real_,
      r2sig        = "",
      slope        = NA_real_,
      intercept    = NA_real_,
      rmse         = NA_real_,
      mae          = NA_real_,
      bias         = NA_real_,
      bias_pct     = NA_real_,
      spearman     = NA_real_,
      spearman_p   = NA_real_,
      spearman_sig = ""
    ))
  }
  
  # Ajuste lineal
  fit <- lm(value ~ local, data = d)
  a <- unname(coef(fit)[1])
  b <- unname(coef(fit)[2])
  
  # Errores
  err <- d$value - d$local
  
  # BIAS (%) = mean((mod-ref)/ref) * 100, evitando ref=0
  rel_bias_vec <- err / d$local
  rel_bias_vec[d$local == 0] <- NA_real_
  bias_pct_val <- mean(rel_bias_vec, na.rm = TRUE) * 100
  
  # Spearman (rho + p) con control robusto
  sp <- tryCatch({
    ok <- is.finite(d$local) & is.finite(d$value)
    x <- d$local[ok]
    y <- d$value[ok]
    
    if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
      list(rho = NA_real_, p = NA_real_)
    } else {
      ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
      list(rho = unname(ct$estimate), p = unname(ct$p.value))
    }
  }, error = function(e) list(rho = NA_real_, p = NA_real_))
  
  tibble(
    dataset      = ds,
    r2           = summary(fit)$r.squared,
    r2sig        = ifelse(summary(fit)$coefficients[2,4] < 0.05, "*", ""),
    slope        = b,
    intercept    = a,
    rmse         = sqrt(mean(err^2, na.rm = TRUE)),
    mae          = mean(abs(err), na.rm = TRUE),
    bias         = mean(err, na.rm = TRUE),
    bias_pct     = bias_pct_val,
    spearman     = sp$rho,
    spearman_p   = sp$p,
    spearman_sig = ifelse(!is.na(sp$p) && sp$p < 0.05, "*", "")
  )
}) %>%
  mutate(dataset = factor(dataset, levels = c("MRBA60","FireCCI51","MCD64A1","GFED5"))) %>%
  arrange(dataset)

# Tabla final (redondeo SOLO aquí)
tabla_stats_mb <- metrics_mb %>%
  transmute(
    Producto        = as.character(dataset),
    `Ecuación`      = paste0("y = ", round(intercept, 2), " + ", round(slope, 2), "·x"),
    `R²`            = paste0(round(r2, 3), r2sig),
    `RMSE (km²)`    = round(rmse, 3),
    `MAE (km²)`     = round(mae, 3),
    `BIAS (km²)`    = round(bias, 3),
    `BIAS (%)`      = round(bias_pct, 2),
    `Spearman (ρ)`  = paste0(round(spearman, 3), spearman_sig),
    `p Spearman`    = signif(spearman_p, 3)
  )

write.csv(tabla_stats_mb,
          file = file.path(dir_plot,
                           paste0("validacion_MapBiomas_Brasil_stats_",
                                  start_year_common, "_", end_year_common, ".csv")),
          row.names = FALSE)

# ==========================================================
# SCATTER
# ==========================================================
# ==========================================================
# SCATTER
# ==========================================================
# ==========================================================
# SCATTER
# ==========================================================
titulo_mb <- paste0(
  "d) MapBiomas (Brazil; ",
  min(anos_comunes),
  "–",
  max(anos_comunes),
  ")"
)

p_mb <- ggplot(
  df_mb_long,
  aes(
    x     = local,
    y     = value,
    color = dataset,
    fill  = dataset
  )
) +
  geom_abline(
    intercept = 0,
    slope     = 1,
    linewidth = 1.1,
    color     = "black",
    linetype  = "longdash"
  ) +
  geom_point(
    shape = 20,
    size  = 5.5,
    alpha = 0.9
  ) +
  geom_smooth(
    method    = "lm",
    se        = TRUE,
    linewidth = 1.0,
    alpha     = 0.25
  ) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  scale_x_continuous(
    labels = scales::number_format(accuracy = 1)
  ) +
  scale_y_continuous(
    labels = scales::number_format(accuracy = 1)
  ) +
  labs(
    title = titulo_mb,
    x = expression("MapBiomas BA (km"^2~yr^{-1}*")"),
    y = expression("Databases BA (km"^2~yr^{-1}*")")
  ) +
  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(
      face  = "bold",
      size  = 35,
      hjust = 0.5
    ),
    axis.title.x = element_text(
      size   = 32,
      margin = margin(t = 12)
    ),
    axis.title.y = element_text(
      size   = 32,
      margin = margin(r = 12)
    ),
    axis.text = element_text(
      size  = 24,
      color = "gray20"
    ),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(
      fill  = NA,
      color = "gray80"
    ),
    legend.position = "none",
    plot.margin = margin(10, 12, 10, 12)
  ) +
  coord_cartesian(clip = "off")

print(p_mb)

ggsave(
  filename = file.path(
    dir_plot,
    paste0(
      "validacion_MapBiomas_Brasil_interior_",
      min(anos_comunes),
      "_",
      max(anos_comunes),
      "_4prod.pdf"
    )
  ),
  plot   = p_mb,
  width  = 12,
  height = 7,
  dpi    = 300,
  device = "pdf"
)

ggsave(
  filename = file.path(
    dir_plot,
    paste0(
      "validacion_MapBiomas_Brasil_interior_",
      min(anos_comunes),
      "_",
      max(anos_comunes),
      "_4prod.jpeg"
    )
  ),
  plot   = p_mb,
  width  = 12,
  height = 7,
  dpi    = 300,
  device = "jpeg"
)

# ==========================================================
# MENSAJE FINAL DE CONTROL (OPCIONAL)
# ==========================================================
cat("\nTabla de métricas guardada en:\n",
    file.path(dir_plot, paste0("validacion_MapBiomas_Brasil_stats_", start_year_common, "_", end_year_common, ".csv")),
    "\nDiagnóstico de errores guardado en:\n",
    file.path(dir_plot, paste0("diagnostico_err_MapBiomas_Brasil_", start_year_common, "_", end_year_common, ".csv")),
    "\n")

