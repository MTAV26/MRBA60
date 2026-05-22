
# ==========================================================
# Limpieza inicial
# ==========================================================
rm(list = ls())
graphics.off()
gc()

# ==========================================================
# Librerías
# ==========================================================
library(terra)
library(ncdf4)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
library(tibble)
library(scales)
library(raster)

# ==========================================================
# Configuración de rutas y modelo
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

# Directorio final de figuras (paper)
dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure4_R1"
dir.create(dir_plot, showWarnings = FALSE, recursive = TRUE)

# Paleta de colores
col_map <- c(
  "MRBA60" = "blue",
  "FireCCI51" = "#E69F00",
  "MCD64A1"   = "brown4",
  "GFED5"     = "gray30"
)

# ==========================================================
# Helper: construir vector de años mensuales
# ==========================================================
make_years_mensual <- function(start_year, ntime) {
  stopifnot(ntime %% 12 == 0)
  rep(start_year:(start_year + ntime/12 - 1), each = 12)
}

# ==========================================================
# Cargar armonizado (FireCCI60) mensual 0.25º
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
# rm(Fire51);
# gc()
# MCD64A1 CMG (2001–2022) -> submuestreamos 2003–2022
# nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/nc_output/MCD64CMQ_BA_2001_2022.nc"
# nc_mcd <- nc_open(nc_path_mcd)
# BA_MCD64 <- ncvar_get(nc_mcd, "Band1") / 1e4
# BA_MCD64 <- BA_MCD64[,,25:264]          # 2003–2022
# nc_close(nc_mcd)
# BA_MCD64[is.na(BA_MCD64)] <- 0
# 
# years_mensual_mcd64 <- make_years_mensual(2003, dim(BA_MCD64)[3])
# years_unicos_mcd64  <- unique(years_mensual_mcd64)

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

# ==========================================================
# Regiones a analizar (SFDL10 anual por región)
# ==========================================================
regions   <- c("Amazonia", "Siberia", "Sahel")
ruta_SFDL <- "/mnt/disco6tb/SFDL10/SFLD_grid"

# Periodo por región (lo que pediste)
start_year_common <- 2003
end_year_by_region <- c(
  "Amazonia" = 2019,
  "Sahel"    = 2019,
  "Siberia"  = 2022
)

# ==========================================================
# Función auxiliar: sumar mensual -> anual por years_unicos
# ==========================================================
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
# Precomputar anual para los productos globales
# (cada uno con SUS propios years_unicos)
# ==========================================================
BA_fc60_anual  <- sumar_mensual_a_anual(BA_final,     years_mensual_fc60,  years_unicos_fc60)
BA_fc51_anual  <- sumar_mensual_a_anual(BA_FireCCI51, years_mensual_fc51,  years_unicos_fc51)
BA_mcd64_anual <- sumar_mensual_a_anual(BA_MCD64,     years_mensual_mcd64, years_unicos_mcd64)
BA_gfed5_anual <- sumar_mensual_a_anual(BA_GFED5,     years_mensual_gfed5, years_unicos_gfed5)

# ==========================================================
# Acumulador para PLOT GLOBAL (si luego quieres combinar)
# ==========================================================
df_list_global <- list()

# ==========================================================
# Bucle por región: leer SFDL10 anual, recortar y scatter
# ==========================================================
for (region in regions) {
  
  end_year_common <- unname(end_year_by_region[region])
  if (is.na(end_year_common)) {
    warning("No hay end_year definido para: ", region, ". Se omite.")
    next
  }
  
  # ---- 1) Localizar y cargar SFDL10 anual (tifs por año) ----
  patron <- paste0("^BA_Area_025_", region, "_(\\d{4})\\.tif$")
  archivos <- list.files(path = ruta_SFDL, pattern = patron, full.names = TRUE)
  if (length(archivos) == 0) {
    warning("No se encontraron archivos SFDL10 para: ", region, ". Se omite.")
    next
  }
  archivos <- archivos[order(archivos)]
  stk <- stack(archivos)
  
  # Extensión y rejilla de SFDL10 (para recorte)
  e <- extent(stk)
  dim_x <- ncol(stk); dim_y <- nrow(stk)
  
  # Años disponibles en SFDL10
  años_sfdl <- as.numeric(sub(patron, "\\1", basename(archivos)))
  
  # ---- 1b) Convertir stack a array [lon, lat, time] y pasar a km² ----
  nt    <- raster::nlayers(stk)
  dim_x <- raster::ncol(stk)
  dim_y <- raster::nrow(stk)
  
  ar <- array(NA_real_, dim = c(dim_x, dim_y, nt))
  
  for (i in seq_len(nt)) {
    mat <- raster::as.matrix(stk[[i]], na.rm = FALSE)
    
    if (!is.matrix(mat)) {
      mat <- matrix(mat, nrow = dim_y, ncol = dim_x, byrow = TRUE)
    } else {
      if (nrow(mat) != dim_y || ncol(mat) != dim_x) {
        mat <- matrix(as.vector(mat), nrow = dim_y, ncol = dim_x, byrow = TRUE)
      }
    }
    
    dimnames(mat) <- NULL
    
    # Volteo vertical (latitudes creciendo hacia arriba)
    mat_flip <- mat[dim_y:1, , drop = FALSE]
    
    # Transponer a [lon, lat]
    slice <- t(mat_flip)
    dimnames(slice) <- NULL
    
    ar[, , i] <- slice
  }
  
  ar <- ar / 1e6  # m2 -> km² (si ya es km², elimina este /1e6)
  
  # ---- 2) Recorte espacial de los productos globales al dominio de la región ----
  xmin <- e@xmin; xmax <- e@xmax; ymin <- e@ymin; ymax <- e@ymax
  ix_lon <- which(lon_glob >= xmin & lon_glob <= xmax)
  ix_lat <- which(lat_glob >= ymin & lat_glob <= ymax)
  
  fc60_reg  <- BA_fc60_anual[ix_lon, ix_lat, , drop = FALSE]
  fc51_reg  <- BA_fc51_anual[ix_lon, ix_lat, , drop = FALSE]
  mcd64_reg <- BA_mcd64_anual[ix_lon, ix_lat, , drop = FALSE]
  gfed5_reg <- BA_gfed5_anual[ix_lon, ix_lat, , drop = FALSE]
  
  # ---- 3) Años comunes por región (intersección real entre TODOS) ----
  anos_comunes <- Reduce(intersect, list(
    años_sfdl,
    years_unicos_fc60,
    years_unicos_fc51,
    years_unicos_mcd64,
    years_unicos_gfed5
  ))
  
  anos_comunes <- sort(anos_comunes)
  anos_comunes <- anos_comunes[anos_comunes >= start_year_common & anos_comunes <= end_year_common]
  
  if (length(anos_comunes) == 0) {
    warning("No hay años comunes ", start_year_common, "–", end_year_common,
            " para: ", region, ". Se omite.")
    next
  }
  
  # Índices por producto (cada uno con su eje temporal)
  idx_sfdl  <- match(anos_comunes, años_sfdl)
  idx_fc60  <- match(anos_comunes, years_unicos_fc60)
  idx_fc51  <- match(anos_comunes, years_unicos_fc51)
  idx_mcd64 <- match(anos_comunes, years_unicos_mcd64)
  idx_gfed5 <- match(anos_comunes, years_unicos_gfed5)
  
  # ---- 4) Sumas anuales regionales (km²) ----
  sum_local <- sum_fc60 <- sum_fc51 <- sum_mcd64 <- sum_gfed5 <- numeric(length(anos_comunes))
  
  for (k in seq_along(anos_comunes)) {
    il  <- idx_sfdl[k]
    i60 <- idx_fc60[k]
    i51 <- idx_fc51[k]
    imc <- idx_mcd64[k]
    igf <- idx_gfed5[k]
    
    sum_local[k] <- sum(ar[,,il],         na.rm = TRUE)  # SFDL10
    sum_fc60[k]  <- sum(fc60_reg[,,i60],  na.rm = TRUE)  # FireCCI60
    sum_fc51[k]  <- sum(fc51_reg[,,i51],  na.rm = TRUE)  # FireCCI51
    sum_mcd64[k] <- sum(mcd64_reg[,,imc], na.rm = TRUE)  # MCD64A1
    sum_gfed5[k] <- sum(gfed5_reg[,,igf], na.rm = TRUE)  # GFED5
  }
  
  df <- data.frame(
    año   = anos_comunes,
    local = sum_local,
    fc60  = sum_fc60,
    fc51  = sum_fc51,
    mcd64 = sum_mcd64,
    gfed5 = sum_gfed5
  )
  
  df_long <- df %>%
    pivot_longer(cols = c(fc60, fc51, mcd64, gfed5),
                 names_to = "dataset", values_to = "value") %>%
    mutate(
      dataset = factor(dataset,
                       levels = c("fc60", "fc51", "mcd64", "gfed5"),
                       labels = c("MRBA60", "FireCCI51", "MCD64A1", "GFED5")),
      region  = region
    )
  
  # ---- 5) Métricas por dataset (km²) + Spearman (ρ) con significancia 95% ----
  groups <- df_long %>% group_split(dataset, .keep = TRUE)
  
  metrics <- purrr::map_dfr(groups, function(d){
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
    
    fit <- lm(value ~ local, data = d)
    a <- unname(coef(fit)[1])
    b <- unname(coef(fit)[2])
    
    # BIAS (%) = mean( (mod-ref)/ref ) * 100, evitando ref=0
    rel_bias_vec <- (d$value - d$local) / d$local
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
    
    # Errores (para RMSE/MAE/BIAS)
    err <- d$value - d$local
    
    tibble(
      dataset      = ds,
      r2           = summary(fit)$r.squared,
      r2sig        = ifelse(summary(fit)$coefficients[2,4] < 0.05, "*", ""),
      slope        = round(b, 2),
      intercept    = round(a, 2),
      rmse         = round(sqrt(mean(err^2, na.rm = TRUE)), 3),
      mae          = round(mean(abs(err), na.rm = TRUE), 3),
      bias         = round(mean(err, na.rm = TRUE), 3),
      bias_pct     = round(bias_pct_val, 2),
      spearman     = sp$rho,
      spearman_p   = sp$p,
      spearman_sig = ifelse(!is.na(sp$p) && sp$p < 0.05, "*", "")
    )
  }) %>%
    mutate(dataset = factor(dataset, levels = c("MRBA60","FireCCI51","MCD64A1","GFED5"))) %>%
    arrange(dataset)
  
  # ---- 6) Tabla de estadísticas y guardado a CSV (incluye Spearman) ----
  tabla_stats <- metrics %>%
    transmute(
      Producto        = as.character(dataset),
      `Ecuación`      = paste0("y = ", intercept, " + ", slope, "·x"),
      `R²`            = paste0(round(r2, 3), r2sig),
      `RMSE (km²)`    = rmse,
      `MAE (km²)`     = mae,
      `BIAS (km²)`    = bias,
      `BIAS (%)`      = bias_pct,
      `Spearman (ρ)`  = paste0(round(spearman, 3), spearman_sig),
      `p Spearman`    = signif(spearman_p, 3)
    )
  
  nombre_csv <- paste0("validacion_SFDL10_", region, "_", start_year_common, "_", end_year_common, "_stats.csv")
  write.csv(tabla_stats,
            file = file.path(dir_plot, nombre_csv),
            row.names = FALSE)
  
  # ---- 7) Títulos por región (dinámicos) ----
  title_prefix <- c(
    "Amazonia" = "a)",
    "Sahel"    = "b)",
    "Siberia"  = "c)"
  )
  
  plot_title <- paste0(
    ifelse(!is.na(title_prefix[region]), title_prefix[region], ""),
    " FireCCISFDL10 (",
    region,
    "; ",
    start_year_common,
    "–",
    end_year_common,
    ")"
  )
  
  # ---- 8) Scatterplot por región ----
  p <- ggplot(
    df_long,
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
      title = plot_title,
      x = expression("SFDL10 BA (km"^2~yr^{-1}*")"),
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
  
  print(p)
  
  ggsave(
    filename = file.path(
      dir_plot,
      paste0(
        "validacion_SFDL10_",
        region,
        "_",
        start_year_common,
        "_",
        end_year_common,
        "_4prod.pdf"
      )
    ),
    plot   = p,
    width  = 12,
    height = 7,
    dpi    = 300,
    device = "pdf"
  )
  
  ggsave(
    filename = file.path(
      dir_plot,
      paste0(
        "validacion_SFDL10_",
        region,
        "_",
        start_year_common,
        "_",
        end_year_common,
        "_4prod.jpeg"
      )
    ),
    plot   = p,
    width  = 12,
    height = 7,
    dpi    = 300,
    device = "jpeg"
  )
}