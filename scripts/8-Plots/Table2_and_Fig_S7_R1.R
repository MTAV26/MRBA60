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
library(grid)
library(abind)
library(fields)

# ==========================================================
# Parámetros temporales (2003–2024)
# ==========================================================
start_year <- 2003
end_year   <- 2024

n_months_target <- (end_year - start_year + 1) * 12  # 264
n_2003_2018 <- (2018 - start_year + 1) * 12          # 192
n_2019_2024 <- (2024 - 2019 + 1) * 12                # 72

idx_2019_2024_in_2003_2024 <- (n_2003_2018 + 1):n_months_target  # 193:264

# ==========================================================
# Configuración de rutas y modelo
# ==========================================================
Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir      <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv  <- paste0(output_dir, "/csv/")
output_dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/"
output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_RData <- paste0(output_dir, "/RData/")

dirs <- c(
  output_dir,
  output_dir_csv,
  output_dir_plot,
  output_dir_plot_rle,
  output_dir_RData
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# ==========================================================
# Cargar longitud y latitud
# ==========================================================
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")

# ==========================================================
# Cargar MRBA60 completo
# ==========================================================
ruta_RData <- file.path(output_dir_RData, "MRBA60_BA_m2_monthly_2003_2024.RData")
load(ruta_RData)

if (!exists("BA_MRBA60")) {
  stop("No encuentro el objeto BA_MRBA60 en MRBA60_BA_m2_monthly_2003_2024.RData")
}

# ==========================================================
# MRBA60 histórico 2003–2018
# Antes se etiquetaba como MRBA60H/FireCCI60.
# A partir de ahora se mostrará SIEMPRE como MRBA60.
# ==========================================================
BA_final <- BA_MRBA60 / 1e6  # km²

if (is.null(dim(BA_final)) || length(dim(BA_final)) != 3) {
  stop(paste0("BA_final no es 3D. dim = ", paste(dim(BA_final), collapse = " x ")))
}

BA_final[is.na(BA_final)] <- 0

BA_final <- BA_final[, , 1:n_2003_2018, drop = FALSE]  # 2003–2018

cat("MRBA60 histórico (2003–2018): ",
    paste(dim(BA_final), collapse = " x "), "\n")

# ==========================================================
# MRBA60 operativo/Sentinel-3 2019–2024
# Antes se etiquetaba como FireCCIS311.
# A partir de ahora se mostrará SIEMPRE como MRBA60.
# ==========================================================
BA_FireS3_raw <- BA_MRBA60[, , 193:264] / 1e6

if (is.null(dim(BA_FireS3_raw)) || length(dim(BA_FireS3_raw)) != 3) {
  stop(paste0("BA_FireS3_raw no es 3D. dim = ",
              paste(dim(BA_FireS3_raw), collapse = " x ")))
}

BA_FireS3_2019_2024 <- BA_FireS3_raw[, , 1:n_2019_2024, drop = FALSE]
BA_FireS3_2019_2024[is.na(BA_FireS3_2019_2024)] <- 0

cat("MRBA60 operativo/S3 (2019–2024): ",
    paste(dim(BA_FireS3_2019_2024), collapse = " x "), "\n")

# Check espacial
if (
  dim(BA_final)[1] != dim(BA_FireS3_2019_2024)[1] ||
  dim(BA_final)[2] != dim(BA_FireS3_2019_2024)[2]
) {
  stop("Dimensiones espaciales no coinciden entre MRBA60 2003–2018 y MRBA60 2019–2024.")
}

# ==========================================================
# Concatenar MRBA60 2003–2024
# ==========================================================
BA_final <- abind(BA_final, BA_FireS3_2019_2024, along = 3)

cat("MRBA60 completo (2003–2024): ",
    paste(dim(BA_final), collapse = " x "), "\n")

if (dim(BA_final)[3] != n_months_target) {
  stop("BA_final no tiene 264 meses tras la extensión.")
}

# ==========================================================
# Mantener objeto auxiliar 2003–2024 con NA antes de 2019
# Solo para control interno si se necesitara.
# No se usará como producto separado en plots/tablas.
# ==========================================================
BA_FireS3_full_2003_2024 <- array(
  NA_real_,
  dim = c(dim(BA_final)[1], dim(BA_final)[2], n_months_target)
)

BA_FireS3_full_2003_2024[, , idx_2019_2024_in_2003_2024] <- BA_FireS3_2019_2024

cat("Objeto auxiliar MRBA60-S3 full (2003–2024, NA antes de 2019): ",
    paste(dim(BA_FireS3_full_2003_2024), collapse = " x "), "\n")

# ==========================================================
# Cargar FireCCI51 2003–2024
# ==========================================================
load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))

BA_FireCCI51 <- f51 / 1e6
rm(f51)
gc()

if (is.null(dim(BA_FireCCI51)) || length(dim(BA_FireCCI51)) != 3) {
  stop(paste0("BA_FireCCI51 no es 3D. dim = ",
              paste(dim(BA_FireCCI51), collapse = " x ")))
}

n_extra <- n_months_target - dim(BA_FireCCI51)[3]

if (n_extra < 0) {
  stop("FireCCI51 tiene más de 264 meses. Revisar índices.")
}

if (n_extra > 0) {
  BA_FireCCI51_extra <- array(
    NA_real_,
    dim = c(dim(BA_FireCCI51)[1], dim(BA_FireCCI51)[2], n_extra)
  )
  BA_FireCCI51 <- abind(BA_FireCCI51, BA_FireCCI51_extra, along = 3)
}

cat("FireCCI51 (2003–2024, NA si faltan meses finales): ",
    paste(dim(BA_FireCCI51), collapse = " x "), "\n")

if (
  dim(BA_FireCCI51)[1] != dim(BA_final)[1] ||
  dim(BA_FireCCI51)[2] != dim(BA_final)[2]
) {
  stop("Dimensiones espaciales no coinciden entre FireCCI51 y MRBA60.")
}

if (dim(BA_FireCCI51)[3] != n_months_target) {
  stop("FireCCI51 no tiene 264 meses.")
}

# ==========================================================
# Cargar MCD64A1 CMG 2003–2024
# ==========================================================
nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/MCD64CMQ_Monthly_2000-2024.nc"

nc_mcd <- nc_open(nc_path_mcd)
BA_MCD64 <- ncvar_get(nc_mcd, "band_data") / 1e6
nc_close(nc_mcd)

if (is.null(dim(BA_MCD64)) || length(dim(BA_MCD64)) != 3) {
  stop(paste0("BA_MCD64 no es 3D. dim = ",
              paste(dim(BA_MCD64), collapse = " x ")))
}

BA_MCD64 <- BA_MCD64[, , 27:290, drop = FALSE]  # 2003–2024

cat("MCD64A1 (2003–2024): ",
    paste(dim(BA_MCD64), collapse = " x "), "\n")

if (
  dim(BA_MCD64)[1] != dim(BA_final)[1] ||
  dim(BA_MCD64)[2] != dim(BA_final)[2]
) {
  stop("Dimensiones espaciales no coinciden entre MCD64A1 y MRBA60.")
}

if (dim(BA_MCD64)[3] != n_months_target) {
  stop("MCD64A1 no tiene 264 meses.")
}

# ==========================================================
# Cargar GFED5 2003–2024
# ==========================================================
nc_path_gfed <- "/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc"

nc_gfed <- nc_open(nc_path_gfed)
BA_GFED5 <- ncvar_get(nc_gfed, "burned_area") / 1e6
nc_close(nc_gfed)

if (is.null(dim(BA_GFED5)) || length(dim(BA_GFED5)) != 3) {
  stop(paste0("BA_GFED5 no es 3D. dim = ",
              paste(dim(BA_GFED5), collapse = " x ")))
}

if (dim(BA_GFED5)[3] >= n_months_target) {
  BA_GFED5 <- BA_GFED5[, , 1:n_months_target, drop = FALSE]
} else {
  n_extra <- n_months_target - dim(BA_GFED5)[3]
  BA_GFED5_extra <- array(
    NA_real_,
    dim = c(dim(BA_GFED5)[1], dim(BA_GFED5)[2], n_extra)
  )
  BA_GFED5 <- abind(BA_GFED5, BA_GFED5_extra, along = 3)
}

cat("GFED5 (ajustado a 2003–2024): ",
    paste(dim(BA_GFED5), collapse = " x "), "\n")

if (
  dim(BA_GFED5)[1] != dim(BA_final)[1] ||
  dim(BA_GFED5)[2] != dim(BA_final)[2]
) {
  stop("Dimensiones espaciales no coinciden entre GFED5 y MRBA60.")
}

if (dim(BA_GFED5)[3] != n_months_target) {
  stop("GFED5 no tiene 264 meses.")
}

# ==========================================================
# Checks finales
# ==========================================================
if (dim(BA_final)[3] != n_months_target) {
  stop("MRBA60 no tiene 264 meses.")
}

if (dim(BA_FireS3_full_2003_2024)[3] != n_months_target) {
  stop("Objeto auxiliar MRBA60-S3 full no tiene 264 meses.")
}

if (dim(BA_FireCCI51)[3] != n_months_target) {
  stop("FireCCI51 no tiene 264 meses.")
}

if (dim(BA_MCD64)[3] != n_months_target) {
  stop("MCD64A1 no tiene 264 meses.")
}

if (dim(BA_GFED5)[3] != n_months_target) {
  stop("GFED5 no tiene 264 meses.")
}

cat("\nOK: Todos los productos están alineados a 2003–2024 y en la misma grilla.\n")

# ==========================================================
# Años de validación
# ==========================================================
anios <- 2017:2024

# ==========================================================
# Funciones de métricas
# ==========================================================
rmse <- function(observed, predicted) {
  sqrt(mean((predicted - observed)^2, na.rm = TRUE))
}

mae <- function(observed, predicted) {
  mean(abs(predicted - observed), na.rm = TRUE)
}

# ==========================================================
# Productos a comparar
# IMPORTANTE:
# MRBA60 incluye 2003–2018 + 2019–2024.
# No se separa ya en MRBA60H / FireCCIS311 en plots ni tablas.
# ==========================================================
productos_overlay <- list(
  MRBA60    = BA_final,
  FireCCI51 = BA_FireCCI51,
  MCD64A1   = BA_MCD64,
  GFED5     = BA_GFED5
)

levels_prod <- c("MRBA60", "FireCCI51", "MCD64A1", "GFED5")

# ==========================================================
# 1) Agregación de datos
# ==========================================================
lista_pixeles_masked    <- list()
lista_pixeles_annualfit <- list()

min_months <- 1
axis_max   <- 750

for (yr in anios) {
  
  message("Procesando año ", yr)
  
  idx_start <- (yr - 2003) * 12 + 1
  
  # ----------------------------------------------------------
  # Landsat mensual del año
  # ----------------------------------------------------------
  carpeta <- file.path(
    "/mnt/disco6tb/Validacion",
    as.character(yr),
    "12_AreaQuemada_porCelda_Raster"
  )
  
  tifs_mensuales <- list.files(
    carpeta,
    pattern = "\\.tif$",
    full.names = TRUE
  )
  
  if (length(tifs_mensuales) != 12) {
    warning("Faltan datos de Landsat para ", yr,
            " (", length(tifs_mensuales), " meses). Se omite.")
    next
  }
  
  r_stack <- rast(sort(tifs_mensuales))
  
  arr_3d <- aperm(as.array(r_stack), c(2, 1, 3)) / 1e6
  
  landsat <- arr_3d[, rev(seq_along(lat)), ]
  landsat[landsat < 0] <- NA
  
  valid_mask <- !is.na(landsat)
  
  for (nombre in names(productos_overlay)) {
    
    BA_producto <- productos_overlay[[nombre]]
    BA_anual_3d <- BA_producto[, , idx_start:(idx_start + 11), drop = FALSE]
    
    # --------------------------------------------------------
    # A) Ajuste alineado por meses con observación Landsat
    # --------------------------------------------------------
    BA_sum_masked <- apply(BA_anual_3d * valid_mask, c(1, 2), sum, na.rm = TRUE)
    L_total       <- apply(landsat, c(1, 2), sum, na.rm = TRUE)
    
    n_valid <- apply(valid_mask, c(1, 2), sum)
    
    BA_sum_masked[n_valid < min_months] <- NA
    L_total[n_valid < min_months]       <- NA
    
    df_masked <- tibble(
      Anio = yr,
      Producto = nombre,
      Landsat = as.vector(L_total),
      Producto_Estimado = as.vector(BA_sum_masked)
    ) %>%
      dplyr::filter(!is.na(Landsat), !is.na(Producto_Estimado))
    
    lista_pixeles_masked[[paste(yr, nombre)]] <- df_masked
    
    # --------------------------------------------------------
    # B) Ajuste anual completo del producto
    # --------------------------------------------------------
    BA_sum_full <- apply(BA_anual_3d, c(1, 2), sum, na.rm = TRUE)
    
    df_annual <- tibble(
      Anio = yr,
      Producto = nombre,
      Landsat = as.vector(L_total),
      Producto_Estimado = as.vector(BA_sum_full)
    ) %>%
      dplyr::filter(
        !is.na(Landsat),
        Landsat >= 0,
        !is.na(Producto_Estimado)
      )
    
    lista_pixeles_annualfit[[paste(yr, nombre)]] <- df_annual
  }
}

# ==========================================================
# Combinar tablas
# ==========================================================
tabla_pixeles <- bind_rows(lista_pixeles_masked)
tabla_pixeles_annual <- bind_rows(lista_pixeles_annualfit)

tabla_pixeles$Producto <- factor(
  tabla_pixeles$Producto,
  levels = levels_prod
)

tabla_pixeles_annual$Producto <- factor(
  tabla_pixeles_annual$Producto,
  levels = levels_prod
)

# ==========================================================
# 1bis) Suma y media de área quemada por dataset
# ==========================================================
res_suma_media_masked <- tabla_pixeles %>%
  dplyr::group_by(Anio, Producto) %>%
  dplyr::summarise(
    n_pix         = dplyr::n(),
    sum_Landsat   = sum(Landsat, na.rm = TRUE),
    sum_Producto  = sum(Producto_Estimado, na.rm = TRUE),
    mean_Landsat  = mean(Landsat, na.rm = TRUE),
    mean_Producto = mean(Producto_Estimado, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n====================\n",
    "Suma y media (ajuste alineado, BA por píxel)\n",
    "====================\n", sep = "")
print(res_suma_media_masked)

res_suma_media_annual <- tabla_pixeles_annual %>%
  dplyr::group_by(Anio, Producto) %>%
  dplyr::summarise(
    n_pix         = dplyr::n(),
    sum_Landsat   = sum(Landsat, na.rm = TRUE),
    sum_Producto  = sum(Producto_Estimado, na.rm = TRUE),
    mean_Landsat  = mean(Landsat, na.rm = TRUE),
    mean_Producto = mean(Producto_Estimado, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n====================\n",
    "Suma y media (ajuste anual, BA por píxel)\n",
    "====================\n", sep = "")
print(res_suma_media_annual)

# ==========================================================
# 2) Métricas principales
# ==========================================================
estadisticos <- tabla_pixeles %>%
  dplyr::group_by(Anio, Producto) %>%
  dplyr::summarise(
    N = sum(!is.na(Landsat) & !is.na(Producto_Estimado)),
    
    Spearman = cor(
      Landsat,
      Producto_Estimado,
      method = "spearman",
      use = "complete.obs"
    ),
    
    p_spearman = cor.test(
      Landsat,
      Producto_Estimado,
      method = "spearman",
      exact = FALSE
    )$p.value,
    
    R2 = cor(
      Landsat,
      Producto_Estimado,
      use = "complete.obs"
    )^2,
    
    RMSE = rmse(Landsat, Producto_Estimado),
    MAE  = mae(Landsat, Producto_Estimado),
    
    BIAS = 100 * (
      sum(Producto_Estimado - Landsat, na.rm = TRUE) /
        sum(Landsat, na.rm = TRUE)
    ),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_spearman_fdr = p.adjust(p_spearman, method = "fdr"),
    Spearman_sig   = ifelse(p_spearman_fdr < 0.05, "*", "")
  )

# ==========================================================
# 2bis) Regresiones por año y producto
# ==========================================================
regresiones_masked <- tabla_pixeles %>%
  dplyr::group_by(Anio, Producto) %>%
  dplyr::summarise(
    modelo = list(lm(Producto_Estimado ~ Landsat, data = dplyr::cur_data())),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    intercept = purrr::map_dbl(modelo, ~ unname(coef(.x)[1])),
    slope     = purrr::map_dbl(modelo, ~ unname(coef(.x)[2]))
  ) %>%
  dplyr::select(Anio, Producto, intercept, slope)

cat("\n====================\n",
    "Ecuaciones de regresión (ajuste alineado)\n",
    "====================\n", sep = "")
print(regresiones_masked)

regresiones_masked %>%
  dplyr::arrange(Anio, Producto) %>%
  dplyr::mutate(
    ecuacion = sprintf(
      "Año %d - %s: y = %.4f * x + %.4f",
      Anio,
      as.character(Producto),
      slope,
      intercept
    )
  ) %>%
  dplyr::pull(ecuacion) %>%
  purrr::walk(~ cat(.x, "\n"))

regresiones_annual <- tabla_pixeles_annual %>%
  dplyr::group_by(Anio, Producto) %>%
  dplyr::summarise(
    modelo = list(lm(Producto_Estimado ~ Landsat, data = dplyr::cur_data())),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    intercept = purrr::map_dbl(modelo, ~ unname(coef(.x)[1])),
    slope     = purrr::map_dbl(modelo, ~ unname(coef(.x)[2]))
  ) %>%
  dplyr::select(Anio, Producto, intercept, slope)

cat("\n====================\n",
    "Ecuaciones de regresión (ajuste anual)\n",
    "====================\n", sep = "")
print(regresiones_annual)

regresiones_annual %>%
  dplyr::arrange(Anio, Producto) %>%
  dplyr::mutate(
    ecuacion = sprintf(
      "Año %d - %s: y = %.4f * x + %.4f",
      Anio,
      as.character(Producto),
      slope,
      intercept
    )
  ) %>%
  dplyr::pull(ecuacion) %>%
  purrr::walk(~ cat(.x, "\n"))

# ==========================================================
# 3) Exportar CSV con estadísticos, sumas/medias y regresión
# ==========================================================
estadisticos_completos_masked <- estadisticos %>%
  dplyr::left_join(regresiones_masked, by = c("Anio", "Producto")) %>%
  dplyr::left_join(res_suma_media_masked, by = c("Anio", "Producto")) %>%
  dplyr::arrange(Anio, factor(Producto, levels = levels_prod))

dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/"

if (!dir.exists(dir_plot)) {
  dir.create(dir_plot, recursive = TRUE)
}

csv_stats_masked <- file.path(dir_plot, "Table2.csv")

write.csv(
  estadisticos_completos_masked,
  csv_stats_masked,
  row.names = FALSE
)

cat("\nArchivo CSV guardado en:\n", csv_stats_masked, "\n")

# ==========================================================
# 4) Plot por año
# ==========================================================
dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure_S7_R1/"

if (!dir.exists(dir_plot)) {
  dir.create(dir_plot, recursive = TRUE)
}

col_map <- c(
  "MRBA60"    = "blue",
  "FireCCI51" = "#E69F00",
  "MCD64A1"   = "brown4",
  "GFED5"     = "gray30"
)

for (yr in sort(unique(tabla_pixeles$Anio))) {
  
  df_plot <- tabla_pixeles %>%
    dplyr::filter(Anio == yr) %>%
    dplyr::mutate(
      Producto = factor(as.character(Producto), levels = levels_prod)
    )
  
  df_plot_ann <- tabla_pixeles_annual %>%
    dplyr::filter(Anio == yr) %>%
    dplyr::mutate(
      Producto = factor(as.character(Producto), levels = levels_prod)
    )
  
  # ----------------------------------------------------------
  # Límites de ejes
  # ----------------------------------------------------------
  if (is.na(axis_max)) {
    
    M1 <- max(
      df_plot$Landsat,
      df_plot$Producto_Estimado,
      na.rm = TRUE
    )
    
    M2 <- if (nrow(df_plot_ann) > 0) {
      max(
        df_plot_ann$Landsat,
        df_plot_ann$Producto_Estimado,
        na.rm = TRUE
      )
    } else {
      0
    }
    
    ax <- max(M1, M2, na.rm = TRUE)
    lims <- c(0, ax)
    
  } else {
    
    ax <- axis_max
    lims <- c(0, axis_max)
  }
  
  # ----------------------------------------------------------
  # Etiqueta del panel
  # ----------------------------------------------------------
  panel_tag <- c(
    `2017` = "a)",
    `2018` = "b)",
    `2019` = "c)",
    `2020` = "d)",
    `2021` = "e)",
    `2022` = "f)",
    `2023` = "g)",
    `2024` = "h)"
  )[as.character(yr)]
  
  if (is.na(panel_tag)) {
    panel_tag <- ""
  }
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  p_pub <- ggplot(
    df_plot,
    aes(
      x = Landsat,
      y = Producto_Estimado,
      color = Producto
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "22",
      linewidth = 0.6,
      color = "grey30"
    ) +
    geom_point(
      shape = 19,
      alpha = 0.6,
      size = 1,
      stroke = 0.6,
      na.rm = TRUE
    ) +
    geom_smooth(
      aes(fill = Producto),
      method = "lm",
      se = TRUE,
      level = 0.95,
      linewidth = 0.9,
      alpha = 0.22,
      linetype = "solid",
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = col_map,
      breaks = levels_prod,
      labels = levels_prod,
      guide = guide_legend(
        override.aes = list(
          shape = 19,
          fill = "white",
          size = 2.0,
          stroke = 0.9
        ),
        keywidth  = unit(8, "mm"),
        keyheight = unit(5, "mm"),
        byrow = TRUE
      )
    ) +
    scale_fill_manual(
      values = scales::alpha(col_map, 0.75),
      guide = "none"
    ) +
    coord_equal(
      xlim = lims,
      ylim = lims,
      expand = FALSE
    ) +
    labs(
      title = sprintf("%s Validation — %s", panel_tag, yr),
      x = "Landsat Annual BA (km²)",
      y = "Databases Annual BA (km²)",
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        colour = "grey85"
      ),
      axis.text = element_text(
        colour = "grey20",
        size = 12
      ),
      axis.title.x = element_text(
        size = 18
      ),
      axis.title.y = element_text(
        size = 18
      ),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.2
      ),
      legend.box.margin = ggplot2::margin(4, 4, 4, 4),
      legend.margin = ggplot2::margin(3, 3, 3, 3),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.key.size = unit(6, "mm"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      )
    )
  
  print(p_pub)
  
  ggsave(
    filename = file.path(dir_plot, sprintf("Fig_S7_%s.pdf", yr)),
    plot = p_pub,
    width = 7,
    height = 7,
    useDingbats = FALSE
  )
  
  ggsave(
    filename = file.path(dir_plot, sprintf("Fig_S7_%s.png", yr)),
    plot = p_pub,
    width = 7,
    height = 7,
    dpi = 300
  )
}

cat("\nProceso finalizado. En plots y tablas, MRBA60H y FireCCIS311 se muestran como MRBA60.\n")




# ==========================================================
# 4) Plot por año + figura conjunta 3x3
# ==========================================================
dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure_S7_R1/"

if (!dir.exists(dir_plot)) {
  dir.create(dir_plot, recursive = TRUE)
}

col_map <- c(
  "MRBA60"    = "blue",
  "FireCCI51" = "#E69F00",
  "MCD64A1"   = "brown4",
  "GFED5"     = "gray30"
)

# Aquí guardaremos todos los plots individuales para la figura 3x3
plot_list <- list()

for (yr in sort(unique(tabla_pixeles$Anio))) {
  
  df_plot <- tabla_pixeles %>%
    dplyr::filter(Anio == yr) %>%
    dplyr::mutate(
      Producto = factor(as.character(Producto), levels = levels_prod)
    )
  
  df_plot_ann <- tabla_pixeles_annual %>%
    dplyr::filter(Anio == yr) %>%
    dplyr::mutate(
      Producto = factor(as.character(Producto), levels = levels_prod)
    )
  
  # ----------------------------------------------------------
  # Límites de ejes
  # ----------------------------------------------------------
  if (is.na(axis_max)) {
    
    M1 <- max(
      df_plot$Landsat,
      df_plot$Producto_Estimado,
      na.rm = TRUE
    )
    
    M2 <- if (nrow(df_plot_ann) > 0) {
      max(
        df_plot_ann$Landsat,
        df_plot_ann$Producto_Estimado,
        na.rm = TRUE
      )
    } else {
      0
    }
    
    ax <- max(M1, M2, na.rm = TRUE)
    lims <- c(0, ax)
    
  } else {
    
    ax <- axis_max
    lims <- c(0, axis_max)
  }
  
  # ----------------------------------------------------------
  # Etiqueta del panel
  # ----------------------------------------------------------
  panel_tag <- c(
    `2017` = "a)",
    `2018` = "b)",
    `2019` = "c)",
    `2020` = "d)",
    `2021` = "e)",
    `2022` = "f)",
    `2023` = "g)",
    `2024` = "h)"
  )[as.character(yr)]
  
  if (is.na(panel_tag)) {
    panel_tag <- ""
  }
  
  # ----------------------------------------------------------
  # Plot individual
  # ----------------------------------------------------------
  p_pub <- ggplot(
    df_plot,
    aes(
      x = Landsat,
      y = Producto_Estimado,
      color = Producto
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "22",
      linewidth = 0.6,
      color = "grey30"
    ) +
    geom_point(
      shape = 19,
      alpha = 0.6,
      size = 1,
      stroke = 0.6,
      na.rm = TRUE
    ) +
    geom_smooth(
      aes(fill = Producto),
      method = "lm",
      se = TRUE,
      level = 0.95,
      linewidth = 0.9,
      alpha = 0.22,
      linetype = "solid",
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = col_map,
      breaks = levels_prod,
      labels = levels_prod,
      guide = guide_legend(
        override.aes = list(
          shape = 19,
          fill = "white",
          size = 2.0,
          stroke = 0.9
        ),
        keywidth  = unit(8, "mm"),
        keyheight = unit(5, "mm"),
        byrow = TRUE
      )
    ) +
    scale_fill_manual(
      values = scales::alpha(col_map, 0.75),
      guide = "none"
    ) +
    coord_equal(
      xlim = lims,
      ylim = lims,
      expand = FALSE
    ) +
    labs(
      title = sprintf("%s Validation — %s", panel_tag, yr),
      x = "Landsat Annual BA (km²)",
      y = "Databases Annual BA (km²)",
      color = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = ggplot2::margin(b = 6)
      ),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        linewidth = 0.25,
        colour = "grey85"
      ),
      axis.text = element_text(
        colour = "grey20",
        size = 12
      ),
      axis.title.x = element_text(
        size = 18
      ),
      axis.title.y = element_text(
        size = 18
      ),
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.background = element_rect(
        fill = "white",
        colour = "black",
        linewidth = 0.2
      ),
      legend.box.margin = ggplot2::margin(4, 4, 4, 4),
      legend.margin = ggplot2::margin(3, 3, 3, 3),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.key.size = unit(6, "mm"),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.6
      )
    )
  
  print(p_pub)
  
  # Guardar plot en lista para figura conjunta
  plot_list[[as.character(yr)]] <- p_pub
  
  # Guardar plot individual
  ggsave(
    filename = file.path(dir_plot, sprintf("Fig_S7_%s.pdf", yr)),
    plot = p_pub,
    width = 7,
    height = 7,
    useDingbats = FALSE
  )
  
  ggsave(
    filename = file.path(dir_plot, sprintf("Fig_S7_%s.png", yr)),
    plot = p_pub,
    width = 7,
    height = 7,
    dpi = 300
  )
}

# ==========================================================
# 5) Figura conjunta 3x3
# ==========================================================

# Panel vacío para completar la matriz 3x3
empty_panel <- ggplot() +
  theme_void()

plot_list_3x3 <- c(
  plot_list[as.character(2017:2024)],
  list(empty_panel)
)

fig_3x3 <- gridExtra::arrangeGrob(
  grobs = plot_list_3x3,
  ncol = 3,
  nrow = 3
)

# Guardar PDF conjunto
ggsave(
  filename = file.path(dir_plot, "Fig_S7_all_years_3x3.pdf"),
  plot = fig_3x3,
  width = 21,
  height = 21,
  useDingbats = FALSE
)

# Guardar PNG conjunto
ggsave(
  filename = file.path(dir_plot, "Fig_S7_all_years_3x3.png"),
  plot = fig_3x3,
  width = 21,
  height = 21,
  dpi = 300
)

cat("\nFigura conjunta 3x3 guardada en:\n",
    file.path(dir_plot, "Fig_S7_all_years_3x3.pdf"), "\n",
    file.path(dir_plot, "Fig_S7_all_years_3x3.png"), "\n")

