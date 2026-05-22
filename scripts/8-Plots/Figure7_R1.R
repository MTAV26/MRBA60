# ==========================================================
# FIGURE 7 - GLOBAL BA TIME SERIES
# MRBA60, FireCCI51, MCD64A1 and GFED5
#
# Output:
#   - One single figure with 5 panels:
#       a) Annual
#       b) DJF
#       c) MAM
#       d) JJA
#       e) SON
#
#   - Layout: 1 column x 5 rows
#   - Lines: BA by product
#   - Grey bars: MRBA60 - FireCCI51 difference
#
# Axis settings:
#   - Panel a:
#       Left Y axis:  0 to 10
#       Right Y axis: 0, 0.7, 1.4
#
#   - Panels b-e:
#       Left Y axis:  0 to 2.5
#       Right Y axis: 0, 0.3, 0.6
#
# Additional style:
#   - Right Y axis text in blue
#   - Main X and Y axis text in black
#   - Legend at the bottom, after panel e
#   - Legend includes the grey bar: MRBA60 - FireCCI51
#   - Titles use "BA" instead of "burned area"
#   - Axis titles and legend closer to the plotting area
#   - No grey background band behind the bars
# ==========================================================


# ==========================================================
# 0) Clean environment
# ==========================================================

rm(list = ls())
graphics.off()
gc()


# ==========================================================
# 1) Libraries
# ==========================================================

# install.packages(c(
#   "ncdf4", "ggplot2", "dplyr", "tidyr",
#   "patchwork", "scales"
# ))

library(ncdf4)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(scales)


# ==========================================================
# 2) Configuration
# ==========================================================

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir       <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv   <- file.path(output_dir, "csv")
output_dir_plot  <- file.path(output_dir, "plot")
output_dir_RData <- file.path(output_dir, "RData")

out_dir <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure7_R1/"

dir.create(output_dir_csv,  recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_plot, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_RData, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,         recursive = TRUE, showWarnings = FALSE)


# ==========================================================
# 3) Time axis
# ==========================================================

dates_full <- seq.Date(
  from = as.Date("2003-01-01"),
  to   = as.Date("2024-12-01"),
  by   = "month"
)

n_time <- length(dates_full)

if (n_time != 264) {
  stop("The date sequence is not 264 months long.")
}

years_full <- 2003:2024


# ==========================================================
# 4) Read MRBA60
# Expected object:
#   BA_MRBA60 in m2/month
#
# Output:
#   BA_MRBA60_km2 in km2/month
# ==========================================================

ruta_RData_MRBA60 <- file.path(
  output_dir_RData,
  "MRBA60_BA_m2_monthly_2003_2024.RData"
)

if (!file.exists(ruta_RData_MRBA60)) {
  stop("MRBA60 file not found: ", ruta_RData_MRBA60)
}

load(ruta_RData_MRBA60)

if (!exists("BA_MRBA60")) {
  stop("Object BA_MRBA60 not found inside MRBA60_BA_m2_monthly_2003_2024.RData.")
}

if (length(dim(BA_MRBA60)) != 3) {
  stop("BA_MRBA60 must be a 3D array [lon, lat, time].")
}

if (dim(BA_MRBA60)[3] != n_time) {
  stop(
    "MRBA60 temporal dimension is not 264. Current dimension: ",
    paste(dim(BA_MRBA60), collapse = " x ")
  )
}

BA_MRBA60_km2 <- BA_MRBA60 / 1e6
BA_MRBA60_km2[is.na(BA_MRBA60_km2)] <- 0

rm(BA_MRBA60)
gc()

message("MRBA60 dimensions: ", paste(dim(BA_MRBA60_km2), collapse = " x "))


# ==========================================================
# 5) Read FireCCI51
# Expected object:
#   f51 in m2/month
#
# Output:
#   BA_FireCCI51_km2 in km2/month
# ==========================================================

ruta_RData_FireCCI51 <- file.path(
  dir_oss,
  "FireCCI51_2003_2024_0.25degree.RData"
)

if (!file.exists(ruta_RData_FireCCI51)) {
  stop("FireCCI51 file not found: ", ruta_RData_FireCCI51)
}

load(ruta_RData_FireCCI51)

if (!exists("f51")) {
  stop("Object f51 not found inside FireCCI51_2003_2024_0.25degree.RData.")
}

if (length(dim(f51)) != 3) {
  stop("f51 must be a 3D array [lon, lat, time].")
}

if (dim(f51)[3] != n_time) {
  stop(
    "FireCCI51 temporal dimension is not 264. Current dimension: ",
    paste(dim(f51), collapse = " x ")
  )
}

BA_FireCCI51_km2 <- f51 / 1e6
BA_FireCCI51_km2[is.na(BA_FireCCI51_km2)] <- 0

rm(f51)
gc()

message("FireCCI51 dimensions: ", paste(dim(BA_FireCCI51_km2), collapse = " x "))


# ==========================================================
# 6) Read MCD64A1
# NetCDF monthly 2000-2024
#
# Expected variable:
#   band_data in m2/month
#
# Extraction:
#   2003-01 to 2024-12 = indices 27:290
#
# Output:
#   BA_MCD64A1_km2 in km2/month
# ==========================================================

nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/MCD64CMQ_Monthly_2000-2024.nc"

if (!file.exists(nc_path_mcd)) {
  stop("MCD64A1 NetCDF file not found: ", nc_path_mcd)
}

nc_mcd <- ncdf4::nc_open(nc_path_mcd)

if (!("band_data" %in% names(nc_mcd$var))) {
  ncdf4::nc_close(nc_mcd)
  stop("Variable 'band_data' not found in MCD64A1 NetCDF.")
}

BA_MCD64_all <- ncdf4::ncvar_get(nc_mcd, "band_data")
ncdf4::nc_close(nc_mcd)

if (length(dim(BA_MCD64_all)) != 3) {
  stop("MCD64A1 band_data must be a 3D array [lon, lat, time].")
}

if (dim(BA_MCD64_all)[3] < 290) {
  stop(
    "MCD64A1 temporal dimension is shorter than expected. Current dimension: ",
    paste(dim(BA_MCD64_all), collapse = " x ")
  )
}

BA_MCD64A1_km2 <- BA_MCD64_all[, , 27:290, drop = FALSE] / 1e6
BA_MCD64A1_km2[is.na(BA_MCD64A1_km2)] <- 0

rm(BA_MCD64_all)
gc()

if (dim(BA_MCD64A1_km2)[3] != n_time) {
  stop(
    "MCD64A1 extracted temporal dimension is not 264. Current dimension: ",
    paste(dim(BA_MCD64A1_km2), collapse = " x ")
  )
}

message("MCD64A1 dimensions 2003-2024: ", paste(dim(BA_MCD64A1_km2), collapse = " x "))


# ==========================================================
# 7) Read GFED5
# Expected variable:
#   burned_area in m2/month
#
# Output:
#   BA_GFED5_km2 in km2/month
# ==========================================================

nc_path_gfed <- "/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc"

if (!file.exists(nc_path_gfed)) {
  stop("GFED5 NetCDF file not found: ", nc_path_gfed)
}

nc_gfed <- ncdf4::nc_open(nc_path_gfed)

if (!("burned_area" %in% names(nc_gfed$var))) {
  ncdf4::nc_close(nc_gfed)
  stop("Variable 'burned_area' not found in GFED5 NetCDF.")
}

BA_GFED5_km2 <- ncdf4::ncvar_get(nc_gfed, "burned_area") / 1e6
ncdf4::nc_close(nc_gfed)

if (length(dim(BA_GFED5_km2)) != 3) {
  stop("GFED5 burned_area must be a 3D array [lon, lat, time].")
}

if (dim(BA_GFED5_km2)[3] != n_time) {
  stop(
    "GFED5 temporal dimension is not 264. Current dimension: ",
    paste(dim(BA_GFED5_km2), collapse = " x ")
  )
}

BA_GFED5_km2[is.na(BA_GFED5_km2)] <- 0

message("GFED5 dimensions: ", paste(dim(BA_GFED5_km2), collapse = " x "))


# ==========================================================
# 8) Optional spatial dimension check
# ==========================================================

dims_list <- list(
  MRBA60    = dim(BA_MRBA60_km2),
  FireCCI51 = dim(BA_FireCCI51_km2),
  MCD64A1   = dim(BA_MCD64A1_km2),
  GFED5     = dim(BA_GFED5_km2)
)

print(dims_list)

same_xy <- all(
  dims_list$MRBA60[1:2]    == dims_list$FireCCI51[1:2],
  dims_list$MRBA60[1:2]    == dims_list$MCD64A1[1:2],
  dims_list$MRBA60[1:2]    == dims_list$GFED5[1:2]
)

if (!same_xy) {
  warning(
    "Spatial dimensions are not identical across products. ",
    "This script only computes global sums, so it can continue, ",
    "but check the grids before doing pixel-level comparisons."
  )
}


# ==========================================================
# 9) Global monthly sums
# ==========================================================

message("Computing global monthly BA totals...")

ba_total_MRBA60 <- sapply(seq_len(n_time), function(m) {
  sum(BA_MRBA60_km2[, , m], na.rm = TRUE)
})

ba_total_FireCCI51 <- sapply(seq_len(n_time), function(m) {
  sum(BA_FireCCI51_km2[, , m], na.rm = TRUE)
})

ba_total_MCD64A1 <- sapply(seq_len(n_time), function(m) {
  sum(BA_MCD64A1_km2[, , m], na.rm = TRUE)
})

ba_total_GFED5 <- sapply(seq_len(n_time), function(m) {
  sum(BA_GFED5_km2[, , m], na.rm = TRUE)
})

df_monthly <- data.frame(
  Date      = dates_full,
  Year      = as.integer(format(dates_full, "%Y")),
  Month     = as.integer(format(dates_full, "%m")),
  MRBA60    = ba_total_MRBA60,
  FireCCI51 = ba_total_FireCCI51,
  MCD64A1   = ba_total_MCD64A1,
  GFED5     = ba_total_GFED5,
  check.names = FALSE
)

write.csv(
  df_monthly,
  file = file.path(out_dir, "Figure7_global_monthly_BA_km2.csv"),
  row.names = FALSE
)

# ==========================================================
# 10) Percentage difference between MRBA60 and FireCCI51
#     Periods:
#       - 2019-2024
#       - 2019-2022
#
# Formula:
#   Difference (%) = ((MRBA60 - FireCCI51) / FireCCI51) * 100
# ==========================================================

message("Computing percentage differences between MRBA60 and FireCCI51...")

df_monthly_diff <- data.frame(
  date          = dates_full,
  year          = as.integer(format(dates_full, "%Y")),
  month         = as.integer(format(dates_full, "%m")),
  MRBA60_km2    = ba_total_MRBA60,
  FireCCI51_km2 = ba_total_FireCCI51
) %>%
  dplyr::mutate(
    diff_km2 = MRBA60_km2 - FireCCI51_km2,
    diff_pct = dplyr::if_else(
      FireCCI51_km2 > 0,
      100 * diff_km2 / FireCCI51_km2,
      NA_real_
    )
  )


# ----------------------------------------------------------
# Function to compute total-period percentage difference
# ----------------------------------------------------------

compute_period_diff <- function(df, year_ini, year_end) {
  
  df_period <- df %>%
    dplyr::filter(year >= year_ini, year <= year_end)
  
  total_MRBA60 <- sum(df_period$MRBA60_km2, na.rm = TRUE)
  total_F51    <- sum(df_period$FireCCI51_km2, na.rm = TRUE)
  
  diff_km2 <- total_MRBA60 - total_F51
  
  diff_pct <- ifelse(
    total_F51 > 0,
    100 * diff_km2 / total_F51,
    NA_real_
  )
  
  data.frame(
    period        = paste0(year_ini, "-", year_end),
    MRBA60_km2    = total_MRBA60,
    FireCCI51_km2 = total_F51,
    diff_km2      = diff_km2,
    diff_pct      = diff_pct
  )
}


# ----------------------------------------------------------
# Compute requested periods
# ----------------------------------------------------------

diff_2019_2024 <- compute_period_diff(
  df = df_monthly_diff,
  year_ini = 2019,
  year_end = 2024
)

diff_2019_2022 <- compute_period_diff(
  df = df_monthly_diff,
  year_ini = 2019,
  year_end = 2022
)

df_period_diff <- dplyr::bind_rows(
  diff_2019_2024,
  diff_2019_2022
)


# ----------------------------------------------------------
# Annual differences for diagnostic purposes
# ----------------------------------------------------------

df_annual_diff <- df_monthly_diff %>%
  dplyr::filter(year >= 2019, year <= 2024) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(
    MRBA60_km2    = sum(MRBA60_km2, na.rm = TRUE),
    FireCCI51_km2 = sum(FireCCI51_km2, na.rm = TRUE),
    diff_km2      = MRBA60_km2 - FireCCI51_km2,
    diff_pct      = dplyr::if_else(
      FireCCI51_km2 > 0,
      100 * diff_km2 / FireCCI51_km2,
      NA_real_
    ),
    .groups = "drop"
  )


# ----------------------------------------------------------
# Print results
# ----------------------------------------------------------

message("Period-level percentage differences:")
print(df_period_diff)

message("Annual percentage differences:")
print(df_annual_diff)


# ==========================================================
# 10) Annual aggregation
# ==========================================================

df_annual <- df_monthly %>%
  group_by(Year) %>%
  summarise(
    MRBA60    = sum(MRBA60, na.rm = TRUE),
    FireCCI51 = sum(FireCCI51, na.rm = TRUE),
    MCD64A1   = sum(MCD64A1, na.rm = TRUE),
    GFED5     = sum(GFED5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Period = "Annual",
    Panel  = "a) Annual",
    Unit   = "yr"
  )

write.csv(
  df_annual,
  file = file.path(out_dir, "Figure7_global_annual_BA_km2.csv"),
  row.names = FALSE
)


# ==========================================================
# 11) Seasonal aggregation
# DJF is labelled by the year of January-February.
#
# Example:
#   Dec 2003 + Jan 2004 + Feb 2004 = DJF 2004
# ==========================================================

seasonal_aggregate_df <- function(df_monthly) {
  
  df <- df_monthly %>%
    mutate(
      Season = case_when(
        Month %in% c(12, 1, 2)  ~ "DJF",
        Month %in% c(3, 4, 5)   ~ "MAM",
        Month %in% c(6, 7, 8)   ~ "JJA",
        Month %in% c(9, 10, 11) ~ "SON",
        TRUE ~ NA_character_
      ),
      SeasonYear = ifelse(Season == "DJF" & Month == 12, Year + 1, Year)
    )
  
  df_season <- df %>%
    group_by(SeasonYear, Season) %>%
    summarise(
      n_months  = n(),
      MRBA60    = sum(MRBA60, na.rm = TRUE),
      FireCCI51 = sum(FireCCI51, na.rm = TRUE),
      MCD64A1   = sum(MCD64A1, na.rm = TRUE),
      GFED5     = sum(GFED5, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_months == 3) %>%
    rename(Year = SeasonYear) %>%
    filter(Year >= 2003, Year <= 2024) %>%
    mutate(
      Period = Season,
      Panel = case_when(
        Season == "DJF" ~ "b) DJF",
        Season == "MAM" ~ "c) MAM",
        Season == "JJA" ~ "d) JJA",
        Season == "SON" ~ "e) SON",
        TRUE ~ NA_character_
      ),
      Unit = "season"
    ) %>%
    arrange(Year, factor(Period, levels = c("DJF", "MAM", "JJA", "SON")))
  
  df_season
}

df_seasonal <- seasonal_aggregate_df(df_monthly)

write.csv(
  df_seasonal,
  file = file.path(out_dir, "Figure7_global_seasonal_BA_km2.csv"),
  row.names = FALSE
)

# ==========================================================
# 12) Plot preparation function
# ==========================================================

prepare_panel_data <- function(df_panel) {
  
  required_cols <- c("Year", "MRBA60", "FireCCI51", "MCD64A1", "GFED5")
  
  missing_cols <- setdiff(required_cols, names(df_panel))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in df_panel: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  df_panel <- as.data.frame(df_panel)
  
  df_lines <- df_panel %>%
    dplyr::select(
      Year,
      MRBA60,
      FireCCI51,
      MCD64A1,
      GFED5
    ) %>%
    tidyr::pivot_longer(
      cols = c(MRBA60, FireCCI51, MCD64A1, GFED5),
      names_to = "Product",
      values_to = "BA_km2"
    ) %>%
    dplyr::mutate(
      BA_Mkm2 = BA_km2 / 1e6,
      Product = factor(
        Product,
        levels = c("MRBA60", "FireCCI51", "MCD64A1", "GFED5")
      )
    )
  
  df_diff <- df_panel %>%
    dplyr::transmute(
      Year,
      Diff_Mkm2 = (MRBA60 - FireCCI51) / 1e6
    )
  
  list(
    lines = df_lines,
    diff  = df_diff
  )
}


# ==========================================================
# 13) Plot function
# ==========================================================

make_ba_panel <- function(df_panel,
                          panel_title,
                          y_lab,
                          show_x_lab = FALSE,
                          show_legend = FALSE,
                          y_limit = NULL,
                          y_breaks = NULL,
                          diff_limit = NULL,
                          diff_breaks = NULL,
                          diff_axis_lab = expression(Delta*"BA MRBA60 - FireCCI51 (Mkm"^2*")")) {
  
  dat <- prepare_panel_data(df_panel)
  
  df_lines <- dat$lines
  df_diff  <- dat$diff
  
  # ----------------------------------------------------------
  # Main left Y axis
  # ----------------------------------------------------------
  
  if (is.null(y_limit)) {
    y_max <- max(df_lines$BA_Mkm2, na.rm = TRUE) * 1.12
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  } else {
    y_max <- y_limit
  }
  
  if (is.null(y_breaks)) {
    y_breaks <- pretty(c(0, y_max), n = 5)
  }
  
  # ----------------------------------------------------------
  # Lower band for MRBA60 - FireCCI51 difference
  # ----------------------------------------------------------
  
  band_frac  <- 0.22
  band_min_y <- 0
  band_max_y <- y_max * band_frac
  
  if (is.null(diff_limit)) {
    diff_abs_max <- max(abs(df_diff$Diff_Mkm2), na.rm = TRUE)
    if (!is.finite(diff_abs_max) || diff_abs_max == 0) diff_abs_max <- 1
    diff_limit <- ceiling(diff_abs_max * 10) / 10
  }
  
  diff_min <- 0
  diff_max <- diff_limit
  
  if (is.null(diff_breaks)) {
    diff_breaks <- pretty(c(diff_min, diff_max), n = 4)
  }
  
  map_diff_to_y <- function(x) {
    x <- pmax(pmin(x, diff_max), diff_min)
    band_min_y + (x - diff_min) * (band_max_y - band_min_y) / (diff_max - diff_min)
  }
  
  inv_map_y_to_diff <- function(y) {
    diff_min + (y - band_min_y) * (diff_max - diff_min) / (band_max_y - band_min_y)
  }
  
  df_diff <- df_diff %>%
    dplyr::mutate(
      Diff_Mkm2_plot = pmax(pmin(Diff_Mkm2, diff_max), diff_min),
      y0   = map_diff_to_y(0),
      y1   = map_diff_to_y(Diff_Mkm2_plot),
      xmin = Year - 0.36,
      xmax = Year + 0.36
    )
  
  x_min <- min(df_panel$Year, na.rm = TRUE)
  x_max <- max(df_panel$Year, na.rm = TRUE)
  
  year_breaks <- seq(2003, 2024, by = 2)
  year_breaks <- year_breaks[year_breaks >= x_min & year_breaks <= x_max]
  
  # ----------------------------------------------------------
  # Plot
  # ----------------------------------------------------------
  
  p <- ggplot2::ggplot() +
    
    ggplot2::geom_hline(
      yintercept = map_diff_to_y(0),
      linewidth = 0.35,
      linetype = "dashed",
      colour = "blue"
    ) +
    
    ggplot2::geom_rect(
      data = df_diff,
      ggplot2::aes(
        xmin = xmin,
        xmax = xmax,
        ymin = pmin(y0, y1),
        ymax = pmax(y0, y1),
        fill = "MRBA60 - FireCCI51"
      ),
      colour = "grey45",
      linewidth = 0.15,
      alpha = 0.75
    ) +
    
    ggplot2::geom_line(
      data = df_lines,
      ggplot2::aes(
        x = Year,
        y = BA_Mkm2,
        colour = Product,
        group = Product
      ),
      linewidth = 0.72
    ) +
    
    ggplot2::geom_point(
      data = df_lines,
      ggplot2::aes(
        x = Year,
        y = BA_Mkm2,
        colour = Product
      ),
      size = 1.55
    ) +
    
    ggplot2::scale_colour_manual(
      values = c(
        "MRBA60"    = "#0072B2",
        "FireCCI51" = "#E69F00",
        "MCD64A1"   = "#8B3A3A",
        "GFED5"     = "grey25"
      ),
      name = NULL
    ) +
    
    ggplot2::scale_fill_manual(
      values = c("MRBA60 - FireCCI51" = "grey65"),
      name = NULL
    ) +
    
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        order = 1,
        override.aes = list(
          linewidth = 0.8,
          size = 2.0
        )
      ),
      fill = ggplot2::guide_legend(
        order = 2,
        override.aes = list(
          colour = "grey45",
          alpha = 0.75
        )
      )
    ) +
    
    ggplot2::scale_x_continuous(
      breaks = year_breaks,
      limits = c(x_min - 0.75, x_max + 0.75),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    
    ggplot2::scale_y_continuous(
      name = y_lab,
      limits = c(0, y_max),
      breaks = y_breaks,
      expand = ggplot2::expansion(mult = c(0, 0.03)),
      sec.axis = ggplot2::sec_axis(
        transform = ~ inv_map_y_to_diff(.),
        name = diff_axis_lab,
        breaks = diff_breaks
      )
    ) +
    
    ggplot2::labs(
      title = panel_title,
      x = if (show_x_lab) "Year" else NULL
    ) +
    
    ggplot2::theme_bw(base_size = 11) +
    
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = 12,
        face = "bold",
        hjust = 0,
        margin = ggplot2::margin(b = 4),
        colour = "black"
      ),
      
      axis.title.y.left = ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(r = 3),
        colour = "black"
      ),
      
      axis.title.y.right = ggplot2::element_text(
        size = 8.5,
        margin = ggplot2::margin(l = 3),
        colour = "blue"
      ),
      
      axis.title.x = ggplot2::element_text(
        size = 10,
        margin = ggplot2::margin(t = 3),
        colour = "black"
      ),
      
      axis.text.x = ggplot2::element_text(
        size = 8.6,
        colour = "black",
        angle = 0,
        vjust = 0.5
      ),
      
      axis.text.y.left = ggplot2::element_text(
        size = 8.6,
        colour = "black"
      ),
      
      axis.text.y.right = ggplot2::element_text(
        size = 8.2,
        colour = "blue"
      ),
      
      panel.grid.major = ggplot2::element_line(
        linewidth = 0.25,
        colour = "grey88"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        linewidth = 0.45,
        colour = "grey35"
      ),
      
      legend.position = if (show_legend) "bottom" else "none",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.text = ggplot2::element_text(size = 10, colour = "black"),
      legend.key.width = grid::unit(1.1, "cm"),
      legend.key.height = grid::unit(0.35, "cm"),
      legend.margin = ggplot2::margin(t = -2, r = 0, b = -2, l = 0),
      legend.box.margin = ggplot2::margin(t = -4, r = 0, b = 0, l = 0),
      
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 3, l = 4)
    )
  
  return(p)
}
p_annual <- make_ba_panel(
  df_panel = df_annual,
  panel_title = "a) Global annual BA",
  y_lab = expression("BA (Mkm"^2~yr^{-1}*")"),
  show_x_lab = FALSE,
  show_legend = TRUE,
  y_limit = 10,
  y_breaks = c(0, 2.5, 5.0, 7.5, 10),
  diff_limit = 1.4,
  diff_breaks = c(0, 0.7, 1.4),
  diff_axis_lab = expression(Delta*"BA (Mkm"^2~yr^{-1}*")")
)

p_djf <- make_ba_panel(
  df_panel = df_seasonal %>% dplyr::filter(Period == "DJF"),
  panel_title = "b) Global seasonal BA - DJF",
  y_lab = expression("BA (Mkm"^2~season^{-1}*")"),
  show_x_lab = FALSE,
  show_legend = FALSE,
  y_limit = 2.5,
  y_breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5),
  diff_limit = 0.6,
  diff_breaks = c(0, 0.3, 0.6),
  diff_axis_lab = expression(Delta*"BA (Mkm"^2~season^{-1}*")")
)

p_mam <- make_ba_panel(
  df_panel = df_seasonal %>% dplyr::filter(Period == "MAM"),
  panel_title = "c) Global seasonal BA - MAM",
  y_lab = expression("BA (Mkm"^2~season^{-1}*")"),
  show_x_lab = FALSE,
  show_legend = FALSE,
  y_limit = 2.5,
  y_breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5),
  diff_limit = 0.6,
  diff_breaks = c(0, 0.3, 0.6),
  diff_axis_lab = expression(Delta*"BA (Mkm"^2~season^{-1}*")")
)

p_jja <- make_ba_panel(
  df_panel = df_seasonal %>% dplyr::filter(Period == "JJA"),
  panel_title = "d) Global seasonal BA - JJA",
  y_lab = expression("BA (Mkm"^2~season^{-1}*")"),
  show_x_lab = FALSE,
  show_legend = FALSE,
  y_limit = 2.5,
  y_breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5),
  diff_limit = 0.6,
  diff_breaks = c(0, 0.3, 0.6),
  diff_axis_lab = expression(Delta*"BA (Mkm"^2~season^{-1}*")")
)

p_son <- make_ba_panel(
  df_panel = df_seasonal %>% dplyr::filter(Period == "SON"),
  panel_title = "e) Global seasonal BA - SON",
  y_lab = expression("BA (Mkm"^2~season^{-1}*")"),
  show_x_lab = TRUE,
  show_legend = FALSE,
  y_limit = 2.5,
  y_breaks = c(0, 0.5, 1.0, 1.5, 2.0, 2.5),
  diff_limit = 0.6,
  diff_breaks = c(0, 0.3, 0.6),
  diff_axis_lab = expression(Delta*"BA (Mkm"^2~season^{-1}*")")
)


# ==========================================================
# 15) Combine panels: 1 column x 5 rows
# Legend at the bottom, after panel e
# ==========================================================

fig_5panels <- p_annual / p_djf / p_mam / p_jja / p_son +
  patchwork::plot_layout(
    ncol = 1,
    heights = c(1.15, 1, 1, 1, 1),
    guides = "collect"
  ) &
  ggplot2::theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "horizontal",
    legend.text = ggplot2::element_text(size = 10, colour = "black"),
    legend.margin = ggplot2::margin(t = -3, r = 0, b = 0, l = 0),
    legend.box.margin = ggplot2::margin(t = -5, r = 0, b = 0, l = 0)
  )
# ==========================================================
# 16) Save figure
# ==========================================================

out_pdf <- file.path(out_dir, "Figure7_R1_Annual_Seasonal_5panels.pdf")
out_png <- file.path(out_dir, "Figure7_R1_Annual_Seasonal_5panels.png")

ggsave(
  filename = out_pdf,
  plot = fig_5panels,
  width = 18,
  height = 24,
  units = "cm",
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = fig_5panels,
  width = 18,
  height = 24,
  units = "cm",
  dpi = 400
)

message("Figure saved:")
message(out_pdf)
message(out_png)

