# ===============================
# CONFIGURACIÓN PREVIA
# ===============================
rm(list = ls())
graphics.off()
gc()

library(scales)
library(gridExtra)
library(grid)
library(reshape2)
library(ggplot2)
library(ggtext)
library(ncdf4)
library(sp)
library(fields)
library(maps)
library(RColorBrewer)
library(dplyr)
library(lubridate)
library(MASS)
library(sf)
library(terra)
library(raster)
library(rworldmap)
library(graticule)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(caret)
library(gplots)
library(dendextend)
library(corrplot)
library(randomForest)
library(tidyr)
library(tibble)
library(ggplot2)
library(cowplot)
library(fastshap)


Modelo <- "B1-MRBA60-2003-2024"
dir_oss <- '/mnt/disco6tb/MRBA60/data/A3_ADJ/'
output_dir <- paste0("/mnt/disco6tb/MRBA60-2/results/", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0("/mnt/disco6tb/MRBA60-2/results/D1-Plots/")
output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_RData <- paste0(output_dir, "/RData/")

# Crear directorios si no existen
dirs <- c(output_dir, output_dir_csv, output_dir_plot, output_dir_plot_rle, output_dir_RData)
for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# ============================================================================
# CARGAR LONGITUD Y LATITUD
# ============================================================================
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")


# ============================================================================
# FECHAS: PERÍODO COMPLETO Y PERÍODO COMÚN
# ============================================================================
dates_full   <- seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month")
dates_common <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ind_common   <- which(dates_full %in% dates_common)

anni <- 2000:2024
mesi <- rep(1:12, length(anni))
fechas <- seq(as.Date("2000-01-01"), as.Date("2024-12-01"), by = "month")
inicio <- which(fechas == as.Date("2003-01-01"))
fin <- which(fechas == as.Date("2024-12-01"))




# load(paste0(output_dir_RData, "BA_", Modelo, "global_BA_FireHarmonized_Full_Filtered.RData"))

load(paste0(output_dir_RData, "BA_MRBA60.RData"))

load(file.path(dir_oss, "FireCCIS311_2019_2024_0.25degree.RData"))
BA_FireS3 <- s3 / 1e6
BA_FireS3[BA_FireS3 == 0] <- NA
rm(s3)
gc()

load(file.path(dir_oss, "FireCCI51_2003_2024_0.25degree.RData"))
BA_Fire51_tot <- f51 / 1e6
BA_Fire51_tot[BA_Fire51_tot == 0] <- NA
rm(f51)
gc()


nrows <- dim(BA_Fire51_tot)[1]
ncols <- dim(BA_Fire51_tot)[2]
time_common <- dim(BA_FireS3)[3]  # número de meses para el período común

# Crear arreglo para BA_FireS3_tot basado en el período común (para mantener la misma estructura)
BA_FireS3_tot <- array(NA, dim = c(nrows, ncols, dim(BA_Fire51_tot)[3]))
BA_FireS3_tot[,,ind_common] <- BA_FireS3
rm(FireS3, BA_FireS3)
gc()
# 

# Vector de fechas para el periodo total (2001-2022)
dates_tot <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
ntime_tot <- length(dates_tot)

BA_Fire51_tot2=BA_Fire51_tot[,,193:264]
BA_FireS3_tot2=BA_FireS3_tot[,,193:264]
# dim(BA_FIRE60)
global_BA_FireHarmonized_full_filtered2=BA_FIRE60[,,193:264]
# dim(BA_FireS3_tot)
# # Series totales calculadas (suma de cada capa, con na.rm = TRUE)
# ba_total_Fire51 <- sapply(1:ntime_tot, function(m) sum(BA_Fire51_tot2[,,m], na.rm = TRUE))
# ba_total_FireS3 <- sapply(1:ntime_tot, function(m) sum(BA_FireS3_tot2[,,m], na.rm = TRUE))
# ba_total_harmonized_af <- sapply(1:ntime_tot, function(m) sum(global_BA_FireHarmonized_full_filtered2[,,m], na.rm = TRUE))
# ## --- Compute performance metrics ---
# y_true    <- ba_total_FireS3                # referencia (S3)
# y_pred_f  <- ba_total_Fire51                # FireCCI51
# y_pred_h  <- ba_total_harmonized_af         # Harmonised
# 
# ## FireCCI51 vs Ref
# bias_f <- mean(y_pred_f - y_true, na.rm = TRUE)
# rmse_f <- sqrt(mean((y_pred_f - y_true)^2, na.rm = TRUE))
# r2_f   <- summary(lm(y_pred_f ~ y_true))$r.squared
# cor_f  <- cor(y_true, y_pred_f, method = "spearman", use = "pairwise")
# 
# ## Harmonised vs Ref
# bias_avg_before <- mean(y_pred_h - y_true, na.rm = TRUE)
# rmse_avg_before <- sqrt(mean((y_pred_h - y_true)^2, na.rm = TRUE))
# r2_avg_before   <- summary(lm(y_pred_h ~ y_true))$r.squared
# cor_avg_before  <- cor(y_true, y_pred_h, method = "spearman", use = "pairwise")
# 
# 
# # Un vector de fechas mensuales desde 2019-01-01 hasta 2022-12-01
# dates_common_plot <- seq(
#   from = as.Date("2019-01-01"),
#   to   = as.Date("2024-12-01"),
#   by   = "1 month"
# )
# 
# series_temporales <- data.frame(
#   tiempo = 1:ntime_tot,
#   BA_FireCCI51 = ba_total_Fire51,
#   BA_total_FireCIIS311 = ba_total_FireS3,
#   BA_Harmonized   = ba_total_harmonized_af
# )
# 
# 

## --- Colors ---
col_harmonised <- "#D55E00"
col_s3         <- "#56B4E9"
col_cci51      <- "#E69F00"

# ## --- Convert to Mkm² ---
# ba_s3_M   <- ba_total_FireS3 
# ba_f51_M  <- ba_total_Fire51 
# ba_har_M  <- ba_total_harmonized_af 
# 
# 
# range(ba_total_FireS3, na.rm=TRUE)
# range(ba_total_Fire51, na.rm=TRUE)
# range(ba_total_harmonized_af, na.rm=TRUE)
# 
# ## --- Long data ---
# df_plot <- data.frame(
#   Date        = dates_common_plot,
#   FireCCIS311 = ba_s3_M,
#   FireCCI51   = ba_f51_M,
#   Harmonised  = ba_har_M
# )
# df_long <- melt(df_plot, id.vars = "Date",
#                 variable.name = "Dataset", value.name = "BurnedArea")
# 
# # ## --- Force English month labels ---
# # old_loc <- Sys.getlocale("LC_TIME")
# # on.exit(try(Sys.setlocale("LC_TIME", old_loc), silent = TRUE), add = TRUE)
# try(Sys.setlocale("LC_TIME", "C"), silent = TRUE)
# 
# ## --- Breaks ---
# # y_max    <- max(df_long$BurnedArea, na.rm = TRUE)
# # y_breaks <- seq(0, ceiling((y_max*1.05)/0.2)*0.2, by = 0.2)
# y_max <- max(df_long$BurnedArea, na.rm = TRUE)
# y_brks <- pretty(c(0, y_max), n = 7)   # ~6–7 ticks, adaptativos
# ## --- Compute performance metrics ---
# y_true    <- ba_total_FireS3
# y_pred_f  <- ba_total_Fire51
# y_pred_h  <- ba_total_harmonized_af
# 
# bias_f <- mean(y_pred_f - y_true, na.rm = TRUE)
# rmse_f <- sqrt(mean((y_pred_f - y_true)^2, na.rm = TRUE))
# r2_f   <- summary(lm(y_pred_f ~ y_true))$r.squared
# cor_f  <- cor(y_true, y_pred_f, method = "spearman", use = "pairwise")
# 
# bias_avg_before <- mean(y_pred_h - y_true, na.rm = TRUE)
# rmse_avg_before <- sqrt(mean((y_pred_h - y_true)^2, na.rm = TRUE))
# r2_avg_before   <- summary(lm(y_pred_h ~ y_true))$r.squared
# cor_avg_before  <- cor(y_true, y_pred_h, method = "spearman", use = "pairwise")
# 
# ## --- Core plot with legend inside top ---
# cols <- c("MRBA60S" = col_s3, "FireCCI51" = col_cci51, "Harmonised" = col_harmonised)
# # 
# # p_core <- ggplot(df_long, aes(Date, BurnedArea, color = Dataset)) +
# #   geom_line(linewidth = 1.6) +
# #   scale_color_manual(values = cols,
# #                      breaks = c("FireCCIS311","FireCCI51","Harmonised")) +
# #   scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y",
# #                expand = expansion(mult = c(0.01, 0.02))) +
# #   scale_y_continuous(breaks = y_breaks, limits = c(0, max(y_breaks))) +
# #   labs(
# #     title    = "Global Burned Area Monthly Time Series",
# #     subtitle = "(January 2019 – December 2022)",
# #     x = NULL,
# #     y = expression("BA (km"^2*")"),
# #     color = NULL
# #   ) +
# #   theme_minimal(base_size = 16) +
# #   theme(
# #     plot.title    = element_text(face = "bold", hjust = 0.5, size = 24,
# #                                  margin = ggplot2::margin(b = 5)),
# #     plot.subtitle = element_text(hjust = 0.5, size = 18, margin = ggplot2::margin(b = 15)),
# #     legend.position = c(0.5, 0.92),
# #     legend.justification = c(0.5, 0.5),
# #     legend.direction = "horizontal",
# #     legend.text = element_text(size = 16),
# #     # legend.background = element_rect(fill = alpha("white", 0.85), color = "grey40"),
# #     axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
# #     axis.text.y = element_text(size = 14),
# #     axis.title.y = element_text(size = 16),
# #     panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),
# #     panel.grid.minor.x = element_blank(),
# #     panel.grid.major.y = element_line(color = "grey80", linewidth = 0.6),
# #     panel.grid.minor.y = element_blank(),
# #     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9)
# #   ) +
# #   guides(color = guide_legend(nrow = 1, byrow = TRUE))
# 
# ## --- Core plot (replace your scale_y_* line) ---
# # p_core <- ggplot(df_long, aes(Date, BurnedArea, color = Dataset)) +
# #   geom_line(linewidth = 1.6) +
# #   scale_color_manual(values = cols,
# #                      breaks = c("FireCCIS311","FireCCI51","Harmonised")) +
# #   scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y",
# #                expand = expansion(mult = c(0.01, 0.02))) +
# #   scale_y_continuous(
# #     breaks = y_brks,
# #     labels = scales::label_number(big.mark = ",")  # 10,000; 250,000; …
# #   ) +
# #   coord_cartesian(ylim = c(0, y_max * 1.05)) +     # límites limpios
# #   labs(
# #     title    = "c) Global Burned Area Monthly Time Series",
# #     subtitle = "(January 2019 – December 2022)",
# #     x = NULL,
# #     y = expression("BA (km"^2*")"),
# #     color = NULL
# #   ) +
# #   theme_minimal(base_size = 16) +
# #   theme(
# #     plot.title    = element_text(face = "bold", hjust = 0.5, size = 24,
# #                                  margin = ggplot2::margin(b = 5)),
# #     plot.subtitle = element_text(hjust = 0.5, size = 18,
# #                                  margin = ggplot2::margin(b = 15)),
# #     legend.position = c(0.5, 0.92),
# #     legend.justification = c(0.5, 0.5),
# #     legend.direction = "horizontal",
# #     legend.text = element_text(size = 16),
# #     axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
# #     axis.text.y = element_text(size = 14),
# #     axis.title.y = element_text(size = 16),
# #     panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),
# #     panel.grid.major.y = element_line(color = "grey80", linewidth = 0.6),
# #     panel.grid.minor   = element_blank(),
# #     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9)
# #   ) +
# #   guides(color = guide_legend(nrow = 1, byrow = TRUE))
# p_core <- ggplot(df_long, aes(Date, BurnedArea, color = Dataset)) +
#   geom_line(linewidth = 1.6) +
#   scale_color_manual(values = cols,
#                      breaks = c("MRBA60S","FireCCI51","Harmonised")) +
#   scale_x_date(date_breaks = "6 months", date_labels = "%b-%Y",
#                expand = expansion(mult = c(0.01, 0.02))) +
#   scale_y_continuous(
#     breaks = y_brks,
#     labels = scales::label_number(big.mark = ",")
#   ) +
#   coord_cartesian(ylim = c(0, y_max * 1.05)) +
#   labs(
#     title    = "c) Global Burned Area Monthly Time Series",
#     subtitle = "(January 2019 – December 2022)",
#     x = NULL,
#     y = expression("BA (km"^2*")"),
#     color = NULL
#   ) +
#   theme_minimal(base_size = 16) +
#   theme(
#     plot.title    = element_text(face = "bold", hjust = 0.5, size = 30,  # antes 24
#                                  margin = ggplot2::margin(b = 6)),
#     plot.subtitle = element_text(hjust = 0.5, size = 22,                 # antes 18
#                                  margin = ggplot2::margin(b = 18)),
#     legend.position = c(0.5, 0.92),
#     legend.justification = c(0.5, 0.5),
#     legend.direction = "horizontal",
#     legend.text = element_text(size = 18),                               # más grande
#     axis.text.x = element_text(angle = 0, hjust = 0.5, size = 16),
#     axis.text.y = element_text(size = 16),
#     axis.title.y = element_text(size = 18),
#     panel.grid.major.x = element_line(color = "grey80", linewidth = 0.5),
#     panel.grid.major.y = element_line(color = "grey80", linewidth = 0.6),
#     panel.grid.minor   = element_blank(),
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9)
#   ) +
#   guides(color = guide_legend(nrow = 1, byrow = TRUE))
# 
# ## --- Metrics table ---
# # metrics <- data.frame(
# #   Model         = c("FireCCI51 vs FireCCIS311", "Harmonised vs FireCCIS311"),
# #   "BIAS (km²)" = sprintf("%.0f", c(bias_f, bias_avg_before )),
# #   "RMSE (km²)" = sprintf("%.0f", c(rmse_f, rmse_avg_before )),
# #   "Cor"     = sprintf("%.2f", c(cor_f, cor_avg_before)),
# #   "R²"      = sprintf("%.2f", c(r2_f, r2_avg_before)),
# #   check.names   = FALSE
# # )
# metrics <- data.frame(
#   Model        = c("FireCCI51 vs MRBA60S", "Harmonised vs MRBA60S"),
#   "BIAS (km²)" = scales::comma(c(bias_f,         bias_avg_before), accuracy = 1),
#   "RMSE (km²)" = scales::comma(c(rmse_f,         rmse_avg_before), accuracy = 1),
#   "Cor"        = sprintf("%.2f", c(cor_f,        cor_avg_before)),
#   "R²"         = sprintf("%.2f", c(r2_f,         r2_avg_before)),
#   check.names  = FALSE
# )
# # tbl_theme <- ttheme_default(
# #   core = list(
# #     fg_params = list(cex = 1.0, hjust = 0.5, x = 0.5),
# #     bg_params = list(fill = "white", col = NA)
# #   ),
# #   colhead = list(
# #     fg_params = list(cex = 1.1, fontface = "bold", hjust = 0.5, x = 0.5),
# #     bg_params = list(fill = "white", col = NA)
# #   )
# # )
# tbl_theme <- ttheme_default(
#   core = list(
#     fg_params = list(cex = 1.5, hjust = 0.5, x = 0.5),   # antes 1.0
#     bg_params = list(fill = "white", col = NA)
#   ),
#   colhead = list(
#     fg_params = list(cex = 1.5, fontface = "bold", hjust = 0.5, x = 0.5), # antes 1.1
#     bg_params = list(fill = "white", col = NA)
#   )
# )
# 
# tbl <- tableGrob(metrics, rows = NULL, theme = tbl_theme)
# 
# ## --- Compact table width ---
# tbl$widths <- unit.c(
#   unit(8, "cm"),  # Model
#   unit(4.0, "cm"),  # BIAS
#   unit(3.0, "cm"),  # RMSE
#   unit(2.0, "cm"),  # Cor
#   unit(2.0, "cm")   # R²
# )
# 
# ## --- Combine plot + table ---
# final_plot <- gridExtra::arrangeGrob(
#   grobs = list(p_core, tbl),
#   nrow = 2,
#   heights = c(4, 1),
#   widths = c(1)
# )
# 
# grid::grid.newpage(); grid::grid.draw(final_plot)
# # dev.off()
# output_dir_plot="/mnt/disco6tb/MRBA60/results/D1-Plots/"
# # ## --- Optional save ---
# # ggsave(file.path(output_dir_plot, "XGLOBAL_2019_2022_time_series_common_filtered.jpeg"),
# #        plot = final_plot, width = 15, height = 6, dpi = 220)
# # ## --- Save as PDF ---
# # ggsave(
# #   file.path(output_dir_plot, "Figure2.pdf"),
# #   plot = final_plot,
# #   width = 20, height = 4, units = "in"   # en pulgadas
# # )
# 
# dev.off()
## =======================
## Librerías
## =======================
suppressPackageStartupMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
})

## =======================
## Vectorización de arrays
## =======================
# (asegúrate de tener en memoria BA_FireS3_tot2, BA_Fire51_tot2, global_BA_FireHarmonized_full_filtered2)
vec_S3  <- as.vector(BA_FireS3_tot2)                           # FireCCIS311 (Ref)
vec_F51 <- as.vector(BA_Fire51_tot2)                           # FireCCI51 (Est en panel a)
vec_HAR <- as.vector(global_BA_FireHarmonized_full_filtered2)  # Harmonised (Est en panel b)

## Máscaras por panel
valid1 <- is.finite(vec_S3) & is.finite(vec_F51)
valid2 <- is.finite(vec_S3) & is.finite(vec_HAR)

## Data frames por panel
df_cells_f <- data.frame(Ref = vec_S3[valid1], Est = vec_F51[valid1])  # Panel a
df_cells_h <- data.frame(Ref = vec_S3[valid2], Est = vec_HAR[valid2])  # Panel b

## =======================
## Referencia común para el inset (gris)
## =======================
ref_vals_all <- vec_S3[is.finite(vec_S3) & vec_S3 > 0]

## (Ejemplo de colores si no los tienes definidos)
# col_cci51 <- "#F2C14E"      # amarillo
# col_harmonised <- "#F28C28"  # naranja

## =======================
## Función del scatter con inset debajo de los puntos
## =======================
build_scatter <- function(df, title_txt, point_col, y_label,
                          ref_vals, ref_col = "grey30") {
  lims <- range(c(df$Ref, df$Est), na.rm = TRUE)
  dx <- diff(lims); dy <- dx  # coord_equal ⇒ misma escala en X e Y
  
  ## ---- Inset data (log10): Ref global (gris) + Est del panel (color)
  est_vals <- df$Est
  dens_df <- rbind(
    data.frame(val = ref_vals, set = "Ref"),
    data.frame(val = est_vals[is.finite(est_vals) & est_vals > 0], set = "Est")
  )
  dens_df <- subset(dens_df, is.finite(val) & val > 0)
  
  have_inset <- nrow(dens_df) >= 50 && length(unique(dens_df$val)) > 10
  
  inset_cols  <- c("Ref" = ref_col, "Est" = point_col)
  inset_fills <- c("Ref" = scales::alpha(ref_col, 0.25),
                   "Est" = scales::alpha(point_col, 0.25))
  
  ## ---- Scatter base (sin capas que tapen)
  p <- ggplot(df, aes(x = Ref, y = Est)) +
    coord_cartesian(xlim = lims, ylim = lims, expand = FALSE, clip = "on") +
    coord_equal() +
    labs(
      title = title_txt,
      x = "FireCCIS311 (BA km²)",
      y = y_label
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 19),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_text(size = 15),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey80", linewidth = 0.3),
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  
  ## ---- Inset como FONDO (arriba-izquierda), más ancho y bajito, sin espacio inicial
  if (have_inset) {
    # límites exactos del eje X del inset (mínimo positivo y máximo), sin expand
    xlims_inset <- range(dens_df$val, na.rm = TRUE)
    xlims_inset[1] <- max(xlims_inset[1], .Machine$double.eps)
    
    g_inset <- ggplot(dens_df, aes(x = val, color = set, fill = set)) +
      geom_density(linewidth = 0.6, alpha = 0.25, adjust = 1) +
      scale_x_continuous(
        trans = "log10",
        limits = xlims_inset,
        breaks = scales::log_breaks(n = 5),
        labels = scales::label_number(scale_cut = scales::cut_short_scale()),
        expand = c(0, 0)  # sin espacio extra a la izquierda
      ) +
      scale_y_continuous(labels = function(x) paste0(round(x * 100), "%"),
                         expand = c(0.02, 0)) +
      scale_color_manual(values = inset_cols, guide = "none") +
      scale_fill_manual(values  = inset_fills, guide = "none") +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text  = element_text(size = 7),
        axis.ticks = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2, color = "grey85"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
        plot.background  = element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(0, 0, 0, 0)
      )
    
    inset_grob <- ggplotGrob(g_inset)
    
    # Coordenadas del inset (más pegado, más ancho y menos alto)
    xmin <- lims[1] + 0.01 * dx
    xmax <- lims[1] + 0.46 * dx   # más ancho
    ymin <- lims[1] + 0.78 * dy   # más arriba (reduce altura)
    ymax <- lims[1] + 0.995 * dy  # casi tocando el borde superior
    
    p <- p + annotation_custom(
      grob = inset_grob, xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax
    )
  }
  
  ## ---- Capas encima (puntos primero; líneas encima en negro)
  p +
    geom_point(color = point_col, size = 0.4, alpha = 0.3) +
    geom_smooth(
      method = "lm", se = TRUE,
      color = "black", fill = "black", alpha = 0.10, linewidth = 0.9
    ) +
    geom_abline(
      slope = 1, intercept = 0,
      color = "black", linetype = "dashed", linewidth = 0.9
    )
}

## =======================
## Construcción de paneles
## =======================
p_cell_f <- build_scatter(
  df_cells_f,
  "a) FireCCI51 vs FireCCIS311\n(Cell–month: 2019–2024)",
  point_col = col_cci51,
  y_label   = "FireCCI51 (BA km²)",
  ref_vals  = ref_vals_all
)

p_cell_h <- build_scatter(
  df_cells_h,
  "b) Harmonised vs FireCCIS311\n(Cell–month: 2019–2024)",
  point_col = col_harmonised,
  y_label   = "Harmonised (BA km²)",
  ref_vals  = ref_vals_all
)

scatter_plots <- (p_cell_f + p_cell_h) 
scatter_plots
# plot_annotation(caption = "Inset: gris = FireCCIS311 (Ref, común); color = Est (panel a: FireCCI51; panel b: Harmonised)")
# ## === Guardado coherente (PDF) ===
# output_dir_plot <- "/mnt/disco6tb/Dropbox/UAH/FireCCI60/ATBD/Figure6/"
# if (!dir.exists(output_dir_plot)) dir.create(output_dir_plot, recursive = TRUE)

ggsave(
  filename = file.path(output_dir_plot, "Figure2.pdf"),
  plot = scatter_plots, width = 12, height = 6, device = cairo_pdf
)
# dev.off()
ggsave(
  filename = file.path(output_dir_plot, "Figure2.jpeg"),
  plot = scatter_plots,
  width = 12,
  height = 6,
  dpi = 300,
  device = "jpeg"
)

dev.off()
## ========= utilidades de métricas =========
compute_metrics <- function(df) {
  df <- subset(df, is.finite(Ref) & is.finite(Est))
  n  <- nrow(df)
  if (n < 2) {
    return(data.frame(
      n = n, RMSE = NA_real_, R2 = NA_real_,
      Spearman = NA_real_, p_spear = NA_real_,
      Intercept = NA_real_, Slope = NA_real_
    ))
  }
  rmse <- sqrt(mean((df$Est - df$Ref)^2))
  fit  <- lm(Est ~ Ref, data = df)
  r2   <- summary(fit)$r.squared
  b0   <- unname(coef(fit)[1])
  b1   <- unname(coef(fit)[2])
  ct   <- suppressWarnings(cor.test(df$Ref, df$Est, method = "spearman", exact = FALSE))
  rho  <- unname(ct$estimate)
  pval <- unname(ct$p.value)
  
  data.frame(
    n = n,
    RMSE = rmse,
    R2 = r2,
    Spearman = rho,
    p_spear = pval,
    Intercept = b0,
    Slope = b1
  )
}

pretty_print_metrics <- function(m, titulo = "") {
  m[] <- lapply(m, as.numeric)  # asegurar numérico para formateo
  cat("\n", titulo, "\n", sep = "")
  cat(" n        :", format(m$n, big.mark = ","), "\n", sep = "")
  cat(" RMSE     :", sprintf("%.4f", m$RMSE), "\n", sep = "")
  cat(" R^2      :", sprintf("%.4f", m$R2), "\n", sep = "")
  cat(" Spearman :", sprintf("rho=%.4f, p=%g", m$Spearman, m$p_spear), "\n", sep = "")
  cat(" Recta    :", sprintf("y = %.4f + %.4f·x", m$Intercept, m$Slope), "\n", sep = "")
}

## ========= calcular y ver en consola =========
met_f <- compute_metrics(df_cells_f)
met_h <- compute_metrics(df_cells_h)

pretty_print_metrics(met_f, "a) FireCCI51 vs FireCCIS311")
pretty_print_metrics(met_h, "b) Harmonised vs FireCCIS311")

## (opcional) ver como data.frame/tabla:
resumen <- rbind(
  cbind(Panel = "a) FireCCI51", met_f),
  cbind(Panel = "b) Harmonised", met_h)
)
print(within(resumen, {
  RMSE <- round(RMSE, 4); R2 <- round(R2, 4)
  Spearman <- round(Spearman, 4); p_spear <- signif(p_spear, 3)
  Intercept <- round(Intercept, 4); Slope <- round(Slope, 4)
}), row.names = FALSE)

