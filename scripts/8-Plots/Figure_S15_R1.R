# ============================================================
# ANALISIS DE SENSIBILIDAD MRBA60H vs FireCCI60 vs FireCCI51
# Periodo armonizado: 2003-2018
# 
# MRBA60H   = modelo entrenado con 6 años, 2003-2024
# FireCCI60 = modelo entrenado con 4 años, 2003-2022
# FireCCI51 = producto ESA original, 2003-2024
# ============================================================

rm(list = ls())
graphics.off()
gc()

# ------------------------------------------------------------
# 0) Librerías
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(fields)
  library(RColorBrewer)
  library(viridis)
  library(scales)
})

# ------------------------------------------------------------
# 1) Rutas
# ------------------------------------------------------------
# file_mrba60h <- "/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60H.RData"
# file_firecci60 <- "/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-3/RData/BA_FireCCI60.RData"

file_mrba60h <- "/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_BA_m2_monthly_2003_2024.RData"
file_firecci60 <- "/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-3/RData/FireCCI60_BA_m2_monthly_2003_2024.RData"




file_firecci51 <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/FireCCI51_2003_2024_0.25degree.RData"
dir_lonlat <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
out_dir <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/FigureS15_R1/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "csv"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "plot"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "RData"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 2) Fechas
# ------------------------------------------------------------
dates_2003_2024 <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_2003_2022 <- seq(as.Date("2003-01-01"), as.Date("2022-12-01"), by = "month")

dates_eval <- seq(as.Date("2003-01-01"), as.Date("2018-12-01"), by = "month")

ind_eval_2024 <- which(dates_2003_2024 %in% dates_eval)
ind_eval_2022 <- which(dates_2003_2022 %in% dates_eval)

stopifnot(length(ind_eval_2024) == 192)
stopifnot(length(ind_eval_2022) == 192)

# ------------------------------------------------------------
# 3) Función auxiliar para cargar arrays 3D
# ------------------------------------------------------------
get_first_3d_array <- function(file, preferred_names = NULL) {
  
  e <- new.env()
  load(file, envir = e)
  objs <- ls(e)
  
  if (!is.null(preferred_names)) {
    for (nm in preferred_names) {
      if (exists(nm, envir = e)) {
        x <- get(nm, envir = e)
        if (is.array(x) && length(dim(x)) == 3) {
          message("Objeto usado desde ", basename(file), ": ", nm)
          return(x)
        }
      }
    }
  }
  
  candidates <- objs[sapply(objs, function(nm) {
    x <- get(nm, envir = e)
    is.array(x) && length(dim(x)) == 3
  })]
  
  if (length(candidates) == 0) {
    stop("No se encontró ningún array 3D en: ", file)
  }
  
  nm <- candidates[1]
  message("Objeto usado desde ", basename(file), ": ", nm)
  get(nm, envir = e)
}

# ------------------------------------------------------------
# 4) Cargar productos
# ------------------------------------------------------------

# MRBA60H: producto entrenado con 6 años
BA_MRBA60H <- get_first_3d_array(
  file_mrba60h,
  preferred_names = c(
    "MRBA60H",
    "BA_MRBA60H",
    "global_BA_FireHarmonized_full",
    "global_BA_FireHarmonized_full_filtered_restored",
    "BA_FireHarmonized_full"
  )
)

# FireCCI60: producto entrenado con 4 años
BA_FireCCI60 <- get_first_3d_array(
  file_firecci60,
  preferred_names = c(
    "BA_FireCCI60",
    "FireCCI60",
    "global_BA_FireHarmonized_full",
    "BA_FireHarmonized_full"
  )
)

# FireCCI51: producto ESA original. Normalmente está guardado como f51 en m2.
e51 <- new.env()
load(file_firecci51, envir = e51)

if (exists("f51", envir = e51)) {
  BA_FireCCI51 <- get("f51", envir = e51) / 1e6
  message("Objeto usado desde FireCCI51: f51 dividido por 1e6 para pasar a km2")
} else {
  BA_FireCCI51 <- get_first_3d_array(
    file_firecci51,
    preferred_names = c("BA_FireCCI51", "BA_Fire51_tot", "FireCCI51")
  )
  message("OJO: FireCCI51 no estaba como f51. Revisa si ya está en km2 o si hay que dividir por 1e6.")
}

rm(e51)
gc()
BA_FireCCI60=BA_FireCCI60/1e6
BA_MRBA60H=BA_MRBA60H/1e6
# ------------------------------------------------------------
# 5) Comprobaciones de dimensión
# ------------------------------------------------------------
cat("\nDimensiones:\n")
cat("MRBA60H:   ", paste(dim(BA_MRBA60H), collapse = " x "), "\n")
cat("FireCCI60: ", paste(dim(BA_FireCCI60), collapse = " x "), "\n")
cat("FireCCI51: ", paste(dim(BA_FireCCI51), collapse = " x "), "\n")

stopifnot(length(dim(BA_MRBA60H)) == 3)
stopifnot(length(dim(BA_FireCCI60)) == 3)
stopifnot(length(dim(BA_FireCCI51)) == 3)

stopifnot(dim(BA_MRBA60H)[3] >= max(ind_eval_2024))
stopifnot(dim(BA_FireCCI60)[3] >= max(ind_eval_2022))
stopifnot(dim(BA_FireCCI51)[3] >= max(ind_eval_2024))

stopifnot(all(dim(BA_MRBA60H)[1:2] == dim(BA_FireCCI60)[1:2]))
stopifnot(all(dim(BA_MRBA60H)[1:2] == dim(BA_FireCCI51)[1:2]))

# ------------------------------------------------------------
# 6) Extraer periodo armonizado 2003-2018
# ------------------------------------------------------------
BA_H_eval   <- BA_MRBA60H[, , ind_eval_2024]
BA_60_eval  <- BA_FireCCI60[, , ind_eval_2022]
BA_51_eval  <- BA_FireCCI51[, , ind_eval_2024]

rm(BA_MRBA60H, BA_FireCCI60, BA_FireCCI51)
gc()

# ------------------------------------------------------------
# 7) Funciones de resumen temporal
# ------------------------------------------------------------
sum_monthly <- function(x) {
  sapply(seq_len(dim(x)[3]), function(tt) {
    sum(x[, , tt], na.rm = TRUE)
  })
}

mean_monthly_nonzero <- function(x) {
  sapply(seq_len(dim(x)[3]), function(tt) {
    z <- x[, , tt]
    z <- z[is.finite(z) & z > 0]
    if (length(z) == 0) return(NA_real_)
    mean(z, na.rm = TRUE)
  })
}

n_burned_pixels_monthly <- function(x) {
  sapply(seq_len(dim(x)[3]), function(tt) {
    z <- x[, , tt]
    sum(is.finite(z) & z > 0, na.rm = TRUE)
  })
}

# ------------------------------------------------------------
# 8) Series mensuales globales
# ------------------------------------------------------------
ts_monthly <- data.frame(
  date = dates_eval,
  year = year(dates_eval),
  month = month(dates_eval),
  MRBA60H_6yr = sum_monthly(BA_H_eval),
  FireCCI60_4yr = sum_monthly(BA_60_eval),
  FireCCI51_original = sum_monthly(BA_51_eval),
  mean_nonzero_MRBA60H_6yr = mean_monthly_nonzero(BA_H_eval),
  mean_nonzero_FireCCI60_4yr = mean_monthly_nonzero(BA_60_eval),
  mean_nonzero_FireCCI51_original = mean_monthly_nonzero(BA_51_eval),
  n_pix_MRBA60H_6yr = n_burned_pixels_monthly(BA_H_eval),
  n_pix_FireCCI60_4yr = n_burned_pixels_monthly(BA_60_eval),
  n_pix_FireCCI51_original = n_burned_pixels_monthly(BA_51_eval)
)

ts_monthly <- ts_monthly %>%
  mutate(
    diff_6yr_minus_4yr = MRBA60H_6yr - FireCCI60_4yr,
    rel_6yr_vs_4yr_pct = 100 * (MRBA60H_6yr - FireCCI60_4yr) / FireCCI60_4yr,
    
    diff_6yr_minus_F51 = MRBA60H_6yr - FireCCI51_original,
    diff_4yr_minus_F51 = FireCCI60_4yr - FireCCI51_original,
    
    absdiff_6yr_F51 = abs(diff_6yr_minus_F51),
    absdiff_4yr_F51 = abs(diff_4yr_minus_F51),
    
    closer_to_F51 = case_when(
      absdiff_6yr_F51 < absdiff_4yr_F51 ~ "MRBA60H_6yr closer to FireCCI51",
      absdiff_6yr_F51 > absdiff_4yr_F51 ~ "FireCCI60_4yr closer to FireCCI51",
      TRUE ~ "equal"
    )
  )

write.csv(
  ts_monthly,
  file.path(out_dir, "csv", "monthly_global_sensitivity_2003_2018.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9) Series anuales globales
# ------------------------------------------------------------
ts_annual <- ts_monthly %>%
  group_by(year) %>%
  summarise(
    MRBA60H_6yr = sum(MRBA60H_6yr, na.rm = TRUE),
    FireCCI60_4yr = sum(FireCCI60_4yr, na.rm = TRUE),
    FireCCI51_original = sum(FireCCI51_original, na.rm = TRUE),
    diff_6yr_minus_4yr = MRBA60H_6yr - FireCCI60_4yr,
    rel_6yr_vs_4yr_pct = 100 * (MRBA60H_6yr - FireCCI60_4yr) / FireCCI60_4yr,
    diff_6yr_minus_F51 = MRBA60H_6yr - FireCCI51_original,
    diff_4yr_minus_F51 = FireCCI60_4yr - FireCCI51_original,
    absdiff_6yr_F51 = abs(diff_6yr_minus_F51),
    absdiff_4yr_F51 = abs(diff_4yr_minus_F51),
    .groups = "drop"
  ) %>%
  mutate(
    closer_to_F51 = case_when(
      absdiff_6yr_F51 < absdiff_4yr_F51 ~ "MRBA60H_6yr closer to FireCCI51",
      absdiff_6yr_F51 > absdiff_4yr_F51 ~ "FireCCI60_4yr closer to FireCCI51",
      TRUE ~ "equal"
    )
  )

write.csv(
  ts_annual,
  file.path(out_dir, "csv", "annual_global_sensitivity_2003_2018.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 10) Resumen estadístico global
# ------------------------------------------------------------
summary_global <- data.frame(
  metric = c(
    "Total BA 2003-2018 MRBA60H_6yr",
    "Total BA 2003-2018 FireCCI60_4yr",
    "Total BA 2003-2018 FireCCI51_original",
    "Difference 6yr - 4yr",
    "Relative difference 6yr vs 4yr (%)",
    "Mean monthly difference 6yr - 4yr",
    "Median monthly difference 6yr - 4yr",
    "RMSE monthly 6yr vs 4yr",
    "Correlation monthly 6yr vs 4yr",
    "Months where 6yr > 4yr",
    "Months where 6yr < 4yr",
    "Months where 6yr closer to FireCCI51",
    "Months where 4yr closer to FireCCI51"
  ),
  value = c(
    sum(ts_monthly$MRBA60H_6yr, na.rm = TRUE),
    sum(ts_monthly$FireCCI60_4yr, na.rm = TRUE),
    sum(ts_monthly$FireCCI51_original, na.rm = TRUE),
    sum(ts_monthly$MRBA60H_6yr, na.rm = TRUE) - sum(ts_monthly$FireCCI60_4yr, na.rm = TRUE),
    100 * (sum(ts_monthly$MRBA60H_6yr, na.rm = TRUE) - sum(ts_monthly$FireCCI60_4yr, na.rm = TRUE)) /
      sum(ts_monthly$FireCCI60_4yr, na.rm = TRUE),
    mean(ts_monthly$diff_6yr_minus_4yr, na.rm = TRUE),
    median(ts_monthly$diff_6yr_minus_4yr, na.rm = TRUE),
    sqrt(mean((ts_monthly$MRBA60H_6yr - ts_monthly$FireCCI60_4yr)^2, na.rm = TRUE)),
    cor(ts_monthly$MRBA60H_6yr, ts_monthly$FireCCI60_4yr, use = "complete.obs", method = "spearman"),
    sum(ts_monthly$MRBA60H_6yr > ts_monthly$FireCCI60_4yr, na.rm = TRUE),
    sum(ts_monthly$MRBA60H_6yr < ts_monthly$FireCCI60_4yr, na.rm = TRUE),
    sum(ts_monthly$closer_to_F51 == "MRBA60H_6yr closer to FireCCI51", na.rm = TRUE),
    sum(ts_monthly$closer_to_F51 == "FireCCI60_4yr closer to FireCCI51", na.rm = TRUE)
  )
)

write.csv(
  summary_global,
  file.path(out_dir, "csv", "summary_global_sensitivity_2003_2018.csv"),
  row.names = FALSE
)

print(summary_global)

# ------------------------------------------------------------
# 11) Plot 1: series mensuales
# ------------------------------------------------------------
df_long <- ts_monthly %>%
  dplyr::select(date, MRBA60H_6yr, FireCCI60_4yr, FireCCI51_original) %>%
  tidyr::pivot_longer(
    cols = -date,
    names_to = "product",
    values_to = "BA_km2"
  )

p1 <- ggplot(df_long, aes(x = date, y = BA_km2, colour = product)) +
  geom_line(linewidth = 0.6) +
  theme_bw(base_size = 12) +
  labs(
    title = "Monthly burned area sensitivity, 2003-2018",
    subtitle = "MRBA60H trained with 6 years vs FireCCI60 trained with 4 years vs FireCCI51 original",
    x = NULL,
    y = "Burned area (km²)",
    colour = NULL
  ) +
  scale_y_continuous(labels = scales::comma) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )
p1
ggsave(
  file.path(out_dir, "plot", "P01_monthly_timeseries_2003_2018.png"),
  p1,
  width = 12,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 12) Plot 2: diferencia mensual 6 años - 4 años
# ------------------------------------------------------------
p2 <- ggplot(ts_monthly, aes(x = date, y = diff_6yr_minus_4yr)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Monthly difference between harmonized products, 2003-2018",
    subtitle = "Positive values indicate more burned area in MRBA60H trained with 6 years",
    x = NULL,
    y = "MRBA60H 6yr - FireCCI60 4yr (km²)"
  ) +
  scale_y_continuous(labels = comma) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "plot", "P02_monthly_difference_6yr_minus_4yr_2003_2018.png"),
  p2,
  width = 12,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 13) Plot 3: diferencia relativa mensual
# ------------------------------------------------------------
p3 <- ggplot(ts_monthly, aes(x = date, y = rel_6yr_vs_4yr_pct)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Relative monthly difference, 2003-2018",
    subtitle = "100 × (MRBA60H 6yr - FireCCI60 4yr) / FireCCI60 4yr",
    x = NULL,
    y = "Relative difference (%)"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "plot", "P03_monthly_relative_difference_6yr_vs_4yr_2003_2018.png"),
  p3,
  width = 12,
  height = 5,
  dpi = 300
)

# # ------------------------------------------------------------
# # 14) Plot 4: series anuales
# # ------------------------------------------------------------
# df_annual_long <- ts_annual %>%
#   dplyr::select(year, MRBA60H_6yr, FireCCI60_4yr, FireCCI51_original) %>%
#   tidyr::pivot_longer(
#     cols = -year,
#     names_to = "product",
#     values_to = "BA_km2"
#   )
# 
# p4 <- ggplot(df_annual_long, aes(x = year, y = BA_km2, colour = product)) +
#   geom_line(linewidth = 0.8) +
#   geom_point(size = 2) +
#   theme_bw(base_size = 12) +
#   labs(
#     title = "Annual burned area sensitivity, 2003-2018",
#     x = NULL,
#     y = "Burned area (km²)",
#     colour = NULL
#   ) +
#   scale_x_continuous(breaks = seq(2003, 2018, 1)) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(
#     legend.position = "bottom",
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   file.path(out_dir, "plot", "P04_annual_timeseries_2003_2018.png"),
#   p4,
#   width = 15,
#   height = 4,
#   dpi = 300
# )
# ------------------------------------------------------------
# 14) Plot 4: series anuales mejorado para paper
# ------------------------------------------------------------

df_annual_long <- ts_annual %>%
  dplyr::select(year, MRBA60H_6yr, FireCCI60_4yr, FireCCI51_original) %>%
  tidyr::pivot_longer(
    cols = -year,
    names_to = "product",
    values_to = "BA_km2"
  ) %>%
  dplyr::mutate(
    product = dplyr::case_when(
      product == "MRBA60H_6yr" ~ "MRBA60H, 6-year training",
      product == "FireCCI60_4yr" ~ "FireCCI60, 4-year training",
      product == "FireCCI51_original" ~ "FireCCI51",
      TRUE ~ product
    ),
    BA_Mkm2 = BA_km2 / 1e6
  )

# Orden de la leyenda
df_annual_long$product <- factor(
  df_annual_long$product,
  levels = c(
    "MRBA60H, 6-year training",
    "FireCCI60, 4-year training",
    "FireCCI51"
  )
)

# Colores manuales
cols_products <- c(
  "MRBA60H, 6-year training" = "#D55E00",
  "FireCCI60, 4-year training" = "#0072B2",
  "FireCCI51" = "grey35"
)

# Tipos de línea
lt_products <- c(
  "MRBA60H, 6-year training" = "solid",
  "FireCCI60, 4-year training" = "solid",
  "FireCCI51" = "dashed"
)

p4 <- ggplot2::ggplot(
  df_annual_long,
  ggplot2::aes(
    x = year,
    y = BA_Mkm2,
    colour = product,
    linetype = product,
    group = product
  )
) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::geom_point(size = 2.4, stroke = 0.4) +
  ggplot2::scale_colour_manual(values = cols_products) +
  ggplot2::scale_linetype_manual(values = lt_products) +
  ggplot2::scale_x_continuous(
    breaks = seq(2003, 2018, 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_number(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0.04, 0.08))
  ) +
  ggplot2::labs(
    title = NULL,
    x = NULL,
    y = expression("Annual burned area (10"^6 * " km"^2 * ")"),
    colour = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_classic(base_size = 13) +
  ggplot2::theme(
    axis.title.y = ggplot2::element_text(
      size = 13,
      colour = "black",
      margin = ggplot2::margin(r = 8)
    ),
    axis.text = ggplot2::element_text(size = 11, colour = "black"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
    axis.line = ggplot2::element_line(linewidth = 0.45, colour = "black"),
    axis.ticks = ggplot2::element_line(linewidth = 0.4, colour = "black"),
    axis.ticks.length = grid::unit(0.18, "cm"),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = ggplot2::element_text(size = 11),
    legend.key.width = grid::unit(1.5, "cm"),
    legend.margin = ggplot2::margin(t = 2, r = 0, b = 0, l = 0),
    
    plot.margin = ggplot2::margin(t = 8, r = 10, b = 6, l = 8)
  ) +
  ggplot2::guides(
    colour = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(linewidth = 1.2, size = 2.6)
    ),
    linetype = ggplot2::guide_legend(
      nrow = 1,
      byrow = TRUE
    )
  )

# Mostrar en pantalla
print(p4)

# Guardar PNG
ggplot2::ggsave(
  filename = file.path(out_dir, "plot", "P04_annual_timeseries_2003_2018_paper.png"),
  plot = p4,
  width = 15,
  height = 4.5,
  dpi = 600
)


# Guardar PDF vectorial
ggplot2::ggsave(
  filename = file.path(out_dir, "plot", "P04_annual_timeseries_2003_2018_paper.pdf"),
  plot = p4,
  width = 15,
  height = 4.5
)
# ------------------------------------------------------------
# 15) Plot 5: barras anuales de diferencia
# ------------------------------------------------------------
p5 <- ggplot(ts_annual, aes(x = year, y = diff_6yr_minus_4yr)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col() +
  theme_bw(base_size = 12) +
  labs(
    title = "Annual difference between harmonized products, 2003-2018",
    subtitle = "Positive values indicate more burned area in MRBA60H trained with 6 years",
    x = NULL,
    y = "MRBA60H 6yr - FireCCI60 4yr (km²)"
  ) +
  scale_x_continuous(breaks = seq(2003, 2018, 1)) +
  scale_y_continuous(labels = comma) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "plot", "P05_annual_difference_6yr_minus_4yr_2003_2018.png"),
  p5,
  width = 12,
  height = 5,
  dpi = 300
)
# 
# # ------------------------------------------------------------
# # 16) Plot 6: dispersión mensual 4 años vs 6 años
# # ------------------------------------------------------------
# p6 <- ggplot(ts_monthly, aes(x = FireCCI60_4yr, y = MRBA60H_6yr)) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   geom_point(alpha = 0.75, size = 2) +
#   theme_bw(base_size = 12) +
#   labs(
#     title = "Monthly agreement between harmonized products, 2003-2018",
#     subtitle = "Dashed line = 1:1",
#     x = "FireCCI60 trained with 4 years (km²)",
#     y = "MRBA60H trained with 6 years (km²)"
#   ) +
#   scale_x_continuous(labels = comma) +
#   scale_y_continuous(labels = comma) +
#   theme(
#     plot.title = element_text(face = "bold")
#   )
# 
# ggsave(
#   file.path(out_dir, "plot", "P06_scatter_monthly_6yr_vs_4yr_2003_2018.png"),
#   p6,
#   width = 7,
#   height = 6,
#   dpi = 300
# )
#


# ------------------------------------------------------------
# 16) Plot 6: dispersión mensual MRBA60 (4yr) vs MRBA60 (6yr)
#     Figura mejorada para paper
# ------------------------------------------------------------

# Paquetes necesarios para este bloque
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(grid)
})

# ------------------------------------------------------------
# 16.1) Preparar datos
# ------------------------------------------------------------
df_scatter <- ts_monthly %>%
  dplyr::transmute(
    date = date,
    year = year,
    month = month,
    MRBA60_4yr = FireCCI60_4yr,
    MRBA60_6yr = MRBA60H_6yr
  ) %>%
  dplyr::filter(
    is.finite(MRBA60_4yr),
    is.finite(MRBA60_6yr)
  )

# ------------------------------------------------------------
# 16.2) Métricas de acuerdo
# ------------------------------------------------------------
n_scatter <- nrow(df_scatter)

r_pearson <- cor(
  df_scatter$MRBA60_4yr,
  df_scatter$MRBA60_6yr,
  use = "complete.obs",
  method = "pearson"
)

r_spearman <- cor(
  df_scatter$MRBA60_4yr,
  df_scatter$MRBA60_6yr,
  use = "complete.obs",
  method = "spearman"
)

r2_pearson <- r_pearson^2

rmse_scatter <- sqrt(mean(
  (df_scatter$MRBA60_6yr - df_scatter$MRBA60_4yr)^2,
  na.rm = TRUE
))

mae_scatter <- mean(
  abs(df_scatter$MRBA60_6yr - df_scatter$MRBA60_4yr),
  na.rm = TRUE
)

bias_scatter <- mean(
  df_scatter$MRBA60_6yr - df_scatter$MRBA60_4yr,
  na.rm = TRUE
)

rel_bias_scatter <- 100 * (
  sum(df_scatter$MRBA60_6yr, na.rm = TRUE) -
    sum(df_scatter$MRBA60_4yr, na.rm = TRUE)
) / sum(df_scatter$MRBA60_4yr, na.rm = TRUE)

summary_scatter <- data.frame(
  comparison = "MRBA60 (6yr) vs MRBA60 (4yr)",
  period = "2003-2018",
  n = n_scatter,
  pearson_r = r_pearson,
  spearman_r = r_spearman,
  r2 = r2_pearson,
  rmse_km2 = rmse_scatter,
  mae_km2 = mae_scatter,
  bias_km2 = bias_scatter,
  relative_bias_pct = rel_bias_scatter
)

write.csv(
  summary_scatter,
  file.path(out_dir, "csv", "P06_scatter_monthly_MRBA60_6yr_vs_4yr_2003_2018_metrics.csv"),
  row.names = FALSE
)

print(summary_scatter)

# ------------------------------------------------------------
# 16.3) Límites comunes para eje X/Y
# ------------------------------------------------------------
xy_range <- range(
  c(df_scatter$MRBA60_4yr, df_scatter$MRBA60_6yr),
  na.rm = TRUE
)

xy_pad <- diff(xy_range) * 0.04

xy_limits <- c(
  floor((xy_range[1] - xy_pad) / 50000) * 50000,
  ceiling((xy_range[2] + xy_pad) / 50000) * 50000
)

# Si quieres forzar límites como en tu figura actual, puedes usar esto:
# xy_limits <- c(200000, 950000)

# ------------------------------------------------------------
# 16.4) Texto de métricas dentro del panel
# ------------------------------------------------------------
metrics_label <- paste0(
  "n = ", n_scatter, "\n",
  "R² = ", sprintf("%.3f", r2_pearson), "\n",
  "Spearman = ", sprintf("%.3f", r_spearman), "\n",
  "RMSE = ", scales::comma(round(rmse_scatter, 0)), " km²\n",
  "Bias = ", scales::comma(round(bias_scatter, 0)), " km²",
  " (", sprintf("%.2f", rel_bias_scatter), "%)"
)

# ------------------------------------------------------------
# 16.5) Figura para paper
# ------------------------------------------------------------
p6 <- ggplot2::ggplot(
  df_scatter,
  ggplot2::aes(x = MRBA60_4yr, y = MRBA60_6yr)
) +
  
  # Línea 1:1
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.75,
    colour = "black"
  ) +
  
  # Puntos mensuales
  ggplot2::geom_point(
    shape = 21,
    size = 2.4,
    stroke = 0.35,
    colour = "black",
    fill = "grey25",
    alpha = 0.78
  ) +
  
  # Métricas
  ggplot2::annotate(
    "label",
    x = xy_limits[1] + 0.035 * diff(xy_limits),
    y = xy_limits[2] - 0.035 * diff(xy_limits),
    label = metrics_label,
    hjust = 0,
    vjust = 1,
    size = 3.8,
    label.size = 0.25,
    label.r = grid::unit(0.10, "lines"),
    fill = "white",
    colour = "black"
  ) +
  
  # Escalas iguales
  ggplot2::scale_x_continuous(
    limits = xy_limits,
    labels = scales::comma,
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_y_continuous(
    limits = xy_limits,
    labels = scales::comma,
    breaks = scales::pretty_breaks(n = 5),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::coord_equal() +
  
  # Etiquetas definitivas
  ggplot2::labs(
    title = "Monthly agreement between MRBA60 sensitivity \nruns (2003-2018)",
    # subtitle = "Dashed line indicates the 1:1 relationship",
    x = expression("MRBA60 (4yr) monthly BA (km"^2*")"),
    y = expression("MRBA60 (6yr) monthly BA (km"^2*")")
  ) +
  
  # Tema limpio para paper
  ggplot2::theme_classic(base_size = 14) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = 17,
      face = "bold",
      colour = "black",
      hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(
      size = 13,
      colour = "black",
      hjust = 0,
      margin = ggplot2::margin(t = 3, b = 8)
    ),
    axis.title.x = ggplot2::element_text(
      size = 14,
      colour = "black",
      margin = ggplot2::margin(t = 8)
    ),
    axis.title.y = ggplot2::element_text(
      size = 14,
      colour = "black",
      margin = ggplot2::margin(r = 8)
    ),
    axis.text = ggplot2::element_text(
      size = 12,
      colour = "black"
    ),
    axis.line = ggplot2::element_line(
      linewidth = 0.45,
      colour = "black"
    ),
    axis.ticks = ggplot2::element_line(
      linewidth = 0.4,
      colour = "black"
    ),
    axis.ticks.length = grid::unit(0.18, "cm"),
    panel.grid.major = ggplot2::element_line(
      linewidth = 0.25,
      colour = "grey88"
    ),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(t = 8, r = 10, b = 8, l = 8)
  )

print(p6)

# ------------------------------------------------------------
# 16.6) Guardar figura
# ------------------------------------------------------------

ggplot2::ggsave(
  filename = file.path(out_dir, "plot", "P06_scatter_monthly_MRBA60_6yr_vs_4yr_2003_2018_paper.png"),
  plot = p6,
  width = 12,
  height = 6.5,
  dpi = 600
)

ggplot2::ggsave(
  filename = file.path(out_dir, "plot", "P06_scatter_monthly_MRBA60_6yr_vs_4yr_2003_2018_paper.pdf"),
  plot = p6,
  width = 12,
  height = 6.5,
  device = cairo_pdf
)
# 
# ggplot2::ggsave(
#   filename = file.path(out_dir, "plot", "P06_scatter_monthly_MRBA60_6yr_vs_4yr_2003_2018_paper.tiff"),
#   plot = p6,
#   width = 7.2,
#   height = 6.5,
#   dpi = 600,
#   compression = "lzw"
# )


# ------------------------------------------------------------
# 17) Métricas mensuales por mes calendario
# ------------------------------------------------------------
monthly_clim <- ts_monthly %>%
  group_by(month) %>%
  summarise(
    MRBA60H_6yr = mean(MRBA60H_6yr, na.rm = TRUE),
    FireCCI60_4yr = mean(FireCCI60_4yr, na.rm = TRUE),
    FireCCI51_original = mean(FireCCI51_original, na.rm = TRUE),
    diff_6yr_minus_4yr = mean(diff_6yr_minus_4yr, na.rm = TRUE),
    rel_6yr_vs_4yr_pct = mean(rel_6yr_vs_4yr_pct, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  monthly_clim,
  file.path(out_dir, "csv", "monthly_climatology_sensitivity_2003_2018.csv"),
  row.names = FALSE
)

p7 <- ggplot(monthly_clim, aes(x = month, y = diff_6yr_minus_4yr)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_col() +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = 1:12) +
  scale_y_continuous(labels = comma) +
  labs(
    title = "Mean monthly sensitivity by calendar month, 2003-2018",
    subtitle = "Positive values indicate more burned area in MRBA60H trained with 6 years",
    x = "Month",
    y = "Mean difference 6yr - 4yr (km²)"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(out_dir, "plot", "P07_calendar_month_sensitivity_2003_2018.png"),
  p7,
  width = 9,
  height = 5,
  dpi = 300
)

# ------------------------------------------------------------
# 18) Mapas medios 2003-2018
#     Tratamos NA como 0 para comparar ausencia/presencia de BA.
# ------------------------------------------------------------
calc_mean_map <- function(x) {
  acc <- matrix(0, nrow = dim(x)[1], ncol = dim(x)[2])
  ntime <- dim(x)[3]

  for (tt in seq_len(ntime)) {
    z <- x[, , tt]
    z[!is.finite(z)] <- 0
    acc <- acc + z
  }

  acc / ntime
}

cat("\nCalculando mapas medios 2003-2018...\n")

mean_H  <- calc_mean_map(BA_H_eval)
mean_60 <- calc_mean_map(BA_60_eval)
mean_51 <- calc_mean_map(BA_51_eval)

diff_H_60 <- mean_H - mean_60
diff_H_51 <- mean_H - mean_51
diff_60_51 <- mean_60 - mean_51

rel_H_60 <- 100 * (mean_H - mean_60) / mean_60
rel_H_60[!is.finite(rel_H_60)] <- NA

# ------------------------------------------------------------
# 19) Cargar lon/lat si existen
# ------------------------------------------------------------
lon <- NULL
lat <- NULL

if (file.exists(file.path(dir_lonlat, "longitude.RData"))) {
  e_lon <- new.env()
  load(file.path(dir_lonlat, "longitude.RData"), envir = e_lon)
  lon <- get(ls(e_lon)[1], envir = e_lon)
}

if (file.exists(file.path(dir_lonlat, "latitude.RData"))) {
  e_lat <- new.env()
  load(file.path(dir_lonlat, "latitude.RData"), envir = e_lat)
  lat <- get(ls(e_lat)[1], envir = e_lat)
}

if (is.null(lon)) lon <- seq(-179.875, 179.875, length.out = dim(mean_H)[1])
if (is.null(lat)) lat <- seq(-89.875, 89.875, length.out = dim(mean_H)[2])

# ------------------------------------------------------------
# 20) Función de mapa global simple
# ------------------------------------------------------------
plot_map_png <- function(z, filename, title, zlim = NULL, palette = NULL, legend_lab = "") {

  if (is.null(palette)) {
    palette <- viridis(100)
  }

  png(filename, width = 2400, height = 1200, res = 200)
  par(mar = c(4, 4, 4, 6))

  image.plot(
    lon, lat, z,
    col = palette,
    zlim = zlim,
    xlab = "Longitude",
    ylab = "Latitude",
    main = title,
    legend.lab = legend_lab
  )
  maps::map("world", add = TRUE, col = "grey25", lwd = 0.5)
  # map("world", add = TRUE, col = "grey25", lwd = 0.5)
  box()
  dev.off()
}

# Paleta divergente centrada aproximadamente en 0
div_cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(101)

max_abs_diff <- quantile(abs(diff_H_60), probs = 0.99, na.rm = TRUE)
max_abs_diff <- as.numeric(max_abs_diff)

plot_map_png(
  mean_H,
  file.path(out_dir, "plot", "P08_map_mean_MRBA60H_6yr_2003_2018.png"),
  "Mean monthly BA, MRBA60H trained with 6 years, 2003-2018",
  palette = viridis(100),
  legend_lab = "km²/month"
)

plot_map_png(
  mean_60,
  file.path(out_dir, "plot", "P09_map_mean_FireCCI60_4yr_2003_2018.png"),
  "Mean monthly BA, FireCCI60 trained with 4 years, 2003-2018",
  palette = viridis(100),
  legend_lab = "km²/month"
)

plot_map_png(
  diff_H_60,
  file.path(out_dir, "plot", "P10_map_difference_MRBA60H_6yr_minus_FireCCI60_4yr_2003_2018.png"),
  "Mean monthly BA difference: MRBA60H 6yr - FireCCI60 4yr, 2003-2018",
  zlim = c(-max_abs_diff, max_abs_diff),
  palette = div_cols,
  legend_lab = "km²/month"
)

# ------------------------------------------------------------
# 21) Guardar mapas y objetos principales
# ------------------------------------------------------------
save(
  ts_monthly,
  ts_annual,
  monthly_clim,
  summary_global,
  mean_H,
  mean_60,
  mean_51,
  diff_H_60,
  diff_H_51,
  diff_60_51,
  rel_H_60,
  lon,
  lat,
  file = file.path(out_dir, "RData", "sensitivity_4yr_vs_6yr_2003_2018.RData")
)

# ------------------------------------------------------------
# 22) Limpieza
# ------------------------------------------------------------
rm(BA_H_eval, BA_60_eval, BA_51_eval)
gc()

cat("\n============================================================\n")
cat("ANALISIS FINALIZADO\n")
cat("Resultados guardados en:\n")
cat(out_dir, "\n")
cat("============================================================\n")
