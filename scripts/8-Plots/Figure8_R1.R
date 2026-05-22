##==============================================================================
## FIGURE 9 - GFED REGIONS
## - Regions: GFED basis regions (Giglio 2006)
## - Products: MRBA60, GFED5, MCD64A1
## - Common period: 2003-2024
## - Regional monthly climatology
## - Spatial Spearman correlation by region and month:
##     * MRBA60 vs GFED5
##     * MRBA60 vs MCD64A1
## - Zeros included
## - Final plot:
##     * rho backgrounds
##     * burned area climatology lines
##     * lower rectangle always plotted in front
## - Outputs:
##     * PDF figure
##     * CSV, RDS and Excel tables
##==============================================================================

rm(list = ls())
graphics.off()
gc()

##==============================================================================
## 0. LIBRARIES
##==============================================================================

library(ncdf4)
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

##==============================================================================
## 1. PATHS AND GENERAL SETTINGS
##==============================================================================

Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir       <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv   <- file.path(output_dir, "csv")
output_dir_plot  <- file.path(output_dir, "plot")
output_dir_RData <- file.path(output_dir, "RData")

out_dir_fig <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure8_R1"
out_dir_tab <- file.path(out_dir_fig, "tablas")

dir.create(output_dir_csv,   recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_plot,  recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_RData, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_fig,     recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tab,     recursive = TRUE, showWarnings = FALSE)

file_regions <- "/mnt/disco6tb/FireCCI60/data_005/old/Mask-RegionsGiglio2006/GFED5_Beta_monthly_2002.nc"
file_mrba60  <- file.path(output_dir_RData, "MRBA60_BA_m2_monthly_2003_2024.RData")
file_gfed5   <- "/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc"
file_mcd64   <- "/mnt/disco6tb/MCD64A1_CMG/MCD64CMQ_Monthly_2000-2024.nc"

fig_pdf <- file.path(out_dir_fig, "Figure8_regionesGiglio.pdf")

col_mrba60 <- "blue"
col_gfed5  <- "grey10"
col_mcd64  <- "brown4"

region_names <- c(
  "1"  = "BONA",
  "2"  = "TENA",
  "3"  = "CEAM",
  "4"  = "NHSA",
  "5"  = "SHSA",
  "6"  = "EURO",
  "7"  = "MIDE",
  "8"  = "NHAF",
  "9"  = "SHAF",
  "10" = "BOAS",
  "11" = "CEAS",
  "12" = "SEAS",
  "13" = "EQAS",
  "14" = "AUST"
)

dates_full <- seq.Date(
  from = as.Date("2003-01-01"),
  to   = as.Date("2024-12-01"),
  by   = "month"
)

n_time <- length(dates_full)

if (n_time != 264) {
  stop("The date sequence is not 264 months long.")
}

years  <- as.numeric(format(dates_full, "%Y"))
months <- as.numeric(format(dates_full, "%m"))

##==============================================================================
## 2. FUNCTIONS
##==============================================================================

check_file <- function(path) {
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
}

slice_to_raster <- function(ba_array, t_index, template_raster) {
  
  r <- raster::raster(
    t(ba_array[, , t_index]),
    xmn = raster::xmin(template_raster),
    xmx = raster::xmax(template_raster),
    ymn = raster::ymin(template_raster),
    ymx = raster::ymax(template_raster),
    crs = raster::crs(template_raster)
  )
  
  raster::flip(r, "y")
}

compute_regional_monthly <- function(ba_array, region_raster, years, months) {
  
  out <- vector("list", length(years))
  
  for (i in seq_along(years)) {
    
    r_ba <- raster::raster(
      t(ba_array[, , i]),
      xmn = raster::xmin(region_raster),
      xmx = raster::xmax(region_raster),
      ymn = raster::ymin(region_raster),
      ymx = raster::ymax(region_raster),
      crs = raster::crs(region_raster)
    )
    
    r_ba <- raster::flip(r_ba, "y")
    
    z <- raster::zonal(
      x = r_ba,
      z = region_raster,
      fun = "sum",
      na.rm = TRUE
    )
    
    out[[i]] <- data.frame(
      region_id = z[, 1],
      BA        = z[, 2],
      year      = years[i],
      month     = months[i]
    )
  }
  
  dplyr::bind_rows(out)
}

monthly_mean_raster <- function(ba_array, time_idx, template_raster) {
  
  s <- raster::stack(
    lapply(
      time_idx,
      function(i) slice_to_raster(ba_array, i, template_raster)
    )
  )
  
  raster::calc(s, fun = function(x) mean(x, na.rm = TRUE))
}

safe_spearman <- function(x, y) {
  
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if (n < 3) {
    return(list(rho = NA_real_, p = NA_real_, n = n))
  }
  
  ct <- tryCatch(
    stats::cor.test(x, y, method = "spearman", exact = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(ct)) {
    return(list(rho = NA_real_, p = NA_real_, n = n))
  }
  
  list(
    rho = unname(ct$estimate),
    p   = unname(ct$p.value),
    n   = n
  )
}

sig_stars <- function(p) {
  
  if (is.na(p)) return("")
  if (p <= 0.001) return("***")
  if (p <= 0.01)  return("**")
  if (p <= 0.05)  return("*")
  
  ""
}

compute_spatial_corr_regions <- function(
    ba1,
    ba2,
    region_raster,
    region_ids,
    region_names,
    months_vec,
    min_n = 30
) {
  
  out <- list()
  
  for (m in 1:12) {
    
    idx_m <- which(months_vec == m)
    
    if (length(idx_m) == 0) {
      next
    }
    
    r1_m <- monthly_mean_raster(ba1, idx_m, region_raster)
    r2_m <- monthly_mean_raster(ba2, idx_m, region_raster)
    
    for (rid in region_ids) {
      
      mask_r <- region_raster == rid
      
      v1 <- raster::getValues(
        raster::mask(r1_m, mask_r, maskvalue = 0)
      )
      
      v2 <- raster::getValues(
        raster::mask(r2_m, mask_r, maskvalue = 0)
      )
      
      ok <- is.finite(v1) & is.finite(v2)
      n_ok <- sum(ok)
      
      if (n_ok < min_n) {
        rho  <- NA_real_
        pval <- NA_real_
      } else {
        res  <- safe_spearman(v1, v2)
        rho  <- res$rho
        pval <- res$p
        n_ok <- res$n
      }
      
      out[[length(out) + 1]] <- data.frame(
        region_id   = rid,
        region_name = region_names[as.character(rid)],
        month       = m,
        rho         = rho,
        p_value     = pval,
        n_pixels    = n_ok,
        sig_star    = sig_stars(pval),
        significant = is.finite(pval) & pval <= 0.05
      )
    }
  }
  
  dplyr::bind_rows(out)
}

save_tables <- function(tables_to_save, out_dir) {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (nm in names(tables_to_save)) {
    
    x <- tables_to_save[[nm]]
    
    if (is.null(x)) {
      next
    }
    
    utils::write.csv(
      x,
      file = file.path(out_dir, paste0(nm, ".csv")),
      row.names = FALSE
    )
    
    saveRDS(
      x,
      file = file.path(out_dir, paste0(nm, ".rds"))
    )
  }
  
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    
    wb <- openxlsx::createWorkbook()
    
    for (nm in names(tables_to_save)) {
      
      x <- tables_to_save[[nm]]
      
      if (is.null(x)) {
        next
      }
      
      openxlsx::addWorksheet(wb, nm)
      openxlsx::writeData(wb, sheet = nm, x = x)
    }
    
    openxlsx::saveWorkbook(
      wb,
      file = file.path(out_dir, "Figure9_tables.xlsx"),
      overwrite = TRUE
    )
    
  } else {
    
    message("Package 'openxlsx' is not installed. Excel file was not created.")
  }
}

##==============================================================================
## 3. READ GFED REGIONS
##==============================================================================

check_file(file_regions)

nc_reg <- ncdf4::nc_open(file_regions)

basisregions <- ncdf4::ncvar_get(nc_reg, "basisregions")
lon_reg      <- ncdf4::ncvar_get(nc_reg, "lon")
lat_reg      <- ncdf4::ncvar_get(nc_reg, "lat")

ncdf4::nc_close(nc_reg)

lat_reg <- rev(lat_reg)

library(sp)

r_regions <- raster::raster(
  t(basisregions[, length(lat_reg):1]),
  xmn = min(lon_reg),
  xmx = max(lon_reg),
  ymn = min(lat_reg),
  ymx = max(lat_reg),
  crs = sp::CRS("+proj=longlat +datum=WGS84")
)

r_regions <- raster::flip(r_regions, "y")
r_regions[r_regions[] == 0] <- NA

##==============================================================================
## 4. READ MRBA60
##==============================================================================

check_file(file_mrba60)

load(file_mrba60)

if (!exists("BA_MRBA60")) {
  stop("Object BA_MRBA60 not found inside: ", file_mrba60)
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

ba_mrba60 <- BA_MRBA60 / 1e6
ba_mrba60[is.na(ba_mrba60)] <- 0

rm(BA_MRBA60)
gc()

message("MRBA60 dimensions: ", paste(dim(ba_mrba60), collapse = " x "))

##==============================================================================
## 5. READ GFED5
##==============================================================================

check_file(file_gfed5)

nc_gfed <- ncdf4::nc_open(file_gfed5)

if (!("burned_area" %in% names(nc_gfed$var))) {
  ncdf4::nc_close(nc_gfed)
  stop("Variable 'burned_area' not found in GFED5 NetCDF.")
}

ba_gfed5 <- ncdf4::ncvar_get(nc_gfed, "burned_area") / 1e6

ncdf4::nc_close(nc_gfed)

if (length(dim(ba_gfed5)) != 3) {
  stop("GFED5 burned_area must be a 3D array [lon, lat, time].")
}

if (dim(ba_gfed5)[3] != n_time) {
  stop(
    "GFED5 temporal dimension is not 264. Current dimension: ",
    paste(dim(ba_gfed5), collapse = " x ")
  )
}

ba_gfed5[is.na(ba_gfed5)] <- 0

message("GFED5 dimensions: ", paste(dim(ba_gfed5), collapse = " x "))

##==============================================================================
## 6. READ MCD64A1
##==============================================================================

check_file(file_mcd64)

nc_mcd <- ncdf4::nc_open(file_mcd64)

if (!("band_data" %in% names(nc_mcd$var))) {
  ncdf4::nc_close(nc_mcd)
  stop("Variable 'band_data' not found in MCD64A1 NetCDF.")
}

ba_mcd64_all <- ncdf4::ncvar_get(nc_mcd, "band_data")

ncdf4::nc_close(nc_mcd)

if (length(dim(ba_mcd64_all)) != 3) {
  stop("MCD64A1 band_data must be a 3D array [lon, lat, time].")
}

if (dim(ba_mcd64_all)[3] < 290) {
  stop(
    "MCD64A1 temporal dimension is shorter than expected. Current dimension: ",
    paste(dim(ba_mcd64_all), collapse = " x ")
  )
}

ba_mcd64a1 <- ba_mcd64_all[, , 27:290, drop = FALSE] / 1e6
ba_mcd64a1[is.na(ba_mcd64a1)] <- 0

rm(ba_mcd64_all)
gc()

if (dim(ba_mcd64a1)[3] != n_time) {
  stop(
    "MCD64A1 extracted temporal dimension is not 264. Current dimension: ",
    paste(dim(ba_mcd64a1), collapse = " x ")
  )
}

message("MCD64A1 dimensions 2003-2024: ", paste(dim(ba_mcd64a1), collapse = " x "))

##==============================================================================
## 7. COMMON PERIOD AND GRID ALIGNMENT
##==============================================================================

ntime <- min(
  dim(ba_mrba60)[3],
  dim(ba_gfed5)[3],
  dim(ba_mcd64a1)[3]
)

time_seq <- seq.Date(
  from = as.Date("2003-01-01"),
  by   = "month",
  length.out = ntime
)

years  <- as.numeric(format(time_seq, "%Y"))
months <- as.numeric(format(time_seq, "%m"))

r_template <- raster::raster(
  t(ba_mrba60[, , 1]),
  xmn = raster::xmin(r_regions),
  xmx = raster::xmax(r_regions),
  ymn = raster::ymin(r_regions),
  ymx = raster::ymax(r_regions),
  crs = raster::crs(r_regions)
)

r_template <- raster::flip(r_template, "y")

r_regions_match <- raster::resample(
  x = r_regions,
  y = r_template,
  method = "ngb"
)

region_ids <- sort(stats::na.omit(unique(raster::getValues(r_regions_match))))

##==============================================================================
## 8. REGIONAL MONTHLY SERIES AND MONTHLY CLIMATOLOGY
##==============================================================================

df_mrba60 <- compute_regional_monthly(
  ba_array      = ba_mrba60[, , 1:ntime, drop = FALSE],
  region_raster = r_regions_match,
  years         = years,
  months        = months
)

df_mrba60$product <- "MRBA60"

df_gfed5 <- compute_regional_monthly(
  ba_array      = ba_gfed5[, , 1:ntime, drop = FALSE],
  region_raster = r_regions_match,
  years         = years,
  months        = months
)

df_gfed5$product <- "GFED5"

df_mcd64a1 <- compute_regional_monthly(
  ba_array      = ba_mcd64a1[, , 1:ntime, drop = FALSE],
  region_raster = r_regions_match,
  years         = years,
  months        = months
)

df_mcd64a1$product <- "MCD64A1"

df_all <- dplyr::bind_rows(
  df_mrba60,
  df_gfed5,
  df_mcd64a1
)

df_all <- df_all %>%
  dplyr::mutate(
    region_name = region_names[as.character(region_id)]
  )

clim <- df_all %>%
  dplyr::filter(year >= 2003, year <= 2024) %>%
  dplyr::group_by(product, region_name, month) %>%
  dplyr::summarise(
    BA = mean(BA, na.rm = TRUE),
    .groups = "drop"
  )

##==============================================================================
## 9. OPTIONAL SEASONALITY PLOT
##==============================================================================

p_season <- ggplot(
  clim,
  aes(x = month, y = BA, color = product, shape = product)
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 1.5) +
  facet_wrap(~ region_name, ncol = 2) +
  scale_y_log10() +
  scale_x_continuous(
    breaks = 1:12,
    labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
  ) +
  scale_color_manual(
    values = c(
      "MRBA60" = col_mrba60,
      "GFED5" = col_gfed5,
      "MCD64A1" = col_mcd64
    ),
    name = NULL
  ) +
  labs(
    x = "Month",
    y = expression("Burned area (km"^2*" month"^{-1}*")"),
    shape = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

print(p_season)

##==============================================================================
## 10. SPATIAL SPEARMAN CORRELATION BY REGION AND MONTH
##==============================================================================

df_corr_gfed <- compute_spatial_corr_regions(
  ba1           = ba_mrba60[, , 1:ntime, drop = FALSE],
  ba2           = ba_gfed5[, , 1:ntime, drop = FALSE],
  region_raster = r_regions_match,
  region_ids    = region_ids,
  region_names  = region_names,
  months_vec    = months,
  min_n         = 30
) %>%
  dplyr::mutate(
    comparison = "MRBA60 vs GFED5"
  )

df_corr_mcd <- compute_spatial_corr_regions(
  ba1           = ba_mrba60[, , 1:ntime, drop = FALSE],
  ba2           = ba_mcd64a1[, , 1:ntime, drop = FALSE],
  region_raster = r_regions_match,
  region_ids    = region_ids,
  region_names  = region_names,
  months_vec    = months,
  min_n         = 30
) %>%
  dplyr::mutate(
    comparison = "MRBA60 vs MCD64A1"
  )

df_corr_all <- dplyr::bind_rows(
  df_corr_gfed,
  df_corr_mcd
)

##==============================================================================
## 11. PREPARE BACKGROUND RECTANGLES FOR FINAL FIGURE
##==============================================================================

min_pos <- suppressWarnings(
  min(clim$BA[is.finite(clim$BA) & clim$BA > 0], na.rm = TRUE)
)

if (!is.finite(min_pos)) {
  stop("There is no BA > 0 in 'clim'. Log10 scale cannot be used.")
}

eps <- min_pos

clim_plot <- clim %>%
  dplyr::mutate(
    BA_plot = ifelse(is.finite(BA) & BA > 0, BA, eps)
  )

rho_min <- 0
rho_max <- 1.0

logBA <- log10(clim_plot$BA_plot)

log_min <- min(logBA[is.finite(logBA)], na.rm = TRUE)
log_max <- max(logBA[is.finite(logBA)], na.rm = TRUE)

if (!is.finite(log_min) || !is.finite(log_max)) {
  stop("Invalid log10(BA_plot) range.")
}

if (log_max == log_min) {
  log_min <- log_min - 1
  log_max <- log_max + 1
}

a <- (log_max - log_min) / (rho_max - rho_min)
b <- log_min - a * rho_min

df_bg_gfed <- df_corr_gfed %>%
  dplyr::mutate(
    rho_clip = ifelse(
      is.finite(rho),
      pmin(pmax(rho, rho_min), rho_max),
      NA_real_
    ),
    ymin = eps,
    ymax = ifelse(
      is.finite(rho_clip),
      10^(a * rho_clip + b),
      NA_real_
    ),
    xmin = month - 0.5,
    xmax = month + 0.5,
    bg_legend = "rho (MRBA60 vs GFED5)",
    ymax = ifelse(is.finite(ymax) & ymax >= ymin, ymax, NA_real_)
  ) %>%
  tibble::as_tibble()

df_bg_mcd <- df_corr_mcd %>%
  dplyr::mutate(
    rho_clip = ifelse(
      is.finite(rho),
      pmin(pmax(rho, rho_min), rho_max),
      NA_real_
    ),
    ymin = eps,
    ymax = ifelse(
      is.finite(rho_clip),
      10^(a * rho_clip + b),
      NA_real_
    ),
    xmin = month - 0.5,
    xmax = month + 0.5,
    bg_legend = "rho (MRBA60 vs MCD64A1)",
    ymax = ifelse(is.finite(ymax) & ymax >= ymin, ymax, NA_real_)
  ) %>%
  tibble::as_tibble()

##==============================================================================
## 12. ORDER RECTANGLES
## - The highest rectangle is plotted behind.
## - The lowest rectangle is plotted in front.
##==============================================================================

df_bg_gfed2 <- df_bg_gfed %>%
  dplyr::left_join(
    df_bg_mcd %>%
      dplyr::select(region_name, month, ymax_mcd = ymax),
    by = c("region_name", "month")
  )

df_bg_mcd2 <- df_bg_mcd %>%
  dplyr::left_join(
    df_bg_gfed %>%
      dplyr::select(region_name, month, ymax_gfed = ymax),
    by = c("region_name", "month")
  )

df_gfed_back <- df_bg_gfed2 %>%
  dplyr::filter(
    is.finite(ymax) &
      (is.na(ymax_mcd) | ymax >= ymax_mcd)
  )

df_gfed_front <- df_bg_gfed2 %>%
  dplyr::filter(
    is.finite(ymax) &
      is.finite(ymax_mcd) &
      ymax < ymax_mcd
  )

df_mcd_back <- df_bg_mcd2 %>%
  dplyr::filter(
    is.finite(ymax) &
      is.finite(ymax_gfed) &
      ymax > ymax_gfed
  )

df_mcd_front <- df_bg_mcd2 %>%
  dplyr::filter(
    is.finite(ymax) &
      (is.na(ymax_gfed) | ymax <= ymax_gfed)
  )


##==============================================================================
## 13. FINAL FIGURE - IMPROVED VERSION
## - Larger legend
## - Larger facet titles
## - Region labels with a), b), c), ...
##==============================================================================

month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

##------------------------------------------------------------------------------
## Facet labels: a), b), c), ...
##------------------------------------------------------------------------------

region_order <- sort(unique(clim_plot$region_name))

region_labels <- setNames(
  paste0(letters[seq_along(region_order)], ") ", region_order),
  region_order
)

p <- ggplot() +
  
  ##--------------------------------------------------------------------------
## Background rectangles: spatial correlation
##--------------------------------------------------------------------------
geom_rect(
  data = df_gfed_back,
  aes(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    fill = bg_legend
  ),
  alpha = 0.78,
  color = NA,
  inherit.aes = FALSE
) +
  
  geom_rect(
    data = df_mcd_back,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.55,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = df_mcd_front,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.55,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = df_gfed_front,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.78,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  ##--------------------------------------------------------------------------
## Burned area lines
##--------------------------------------------------------------------------
geom_line(
  data = clim_plot,
  aes(
    x = month,
    y = BA_plot,
    color = product,
    group = product
  ),
  linewidth = 0.82,
  lineend = "round"
) +
  
  geom_point(
    data = clim_plot,
    aes(
      x = month,
      y = BA_plot,
      color = product,
      shape = product
    ),
    size = 2.0,
    stroke = 0.20
  ) +
  
  ##--------------------------------------------------------------------------
## Facets with labels
##--------------------------------------------------------------------------
facet_wrap(
  ~ region_name,
  ncol = 2,
  labeller = as_labeller(region_labels)
) +
  
  ##--------------------------------------------------------------------------
## Axes
##--------------------------------------------------------------------------
scale_y_log10(
  name = expression("Burned area (km"^2*" month"^{-1}*")"),
  breaks = scales::log_breaks(n = 5),
  labels = scales::label_number(
    big.mark = " ",
    decimal.mark = "."
  ),
  sec.axis = sec_axis(
    transform = ~ (log10(.) - b) / a,
    breaks = seq(0, 1.0, by = 0.2),
    labels = function(x) sprintf("%.1f", x),
    name = "Spatial correlation"
  )
) +
  
  scale_x_continuous(
    breaks = 1:12,
    labels = month_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  ##--------------------------------------------------------------------------
## Manual scales
##--------------------------------------------------------------------------
scale_color_manual(
  values = c(
    "MRBA60"  = "#0033FF",
    "GFED5"   = "#202020",
    "MCD64A1" = "#A32020"
  ),
  name = "Burned area"
) +
  
  scale_shape_manual(
    values = c(
      "MRBA60"  = 15,
      "GFED5"   = 16,
      "MCD64A1" = 17
    ),
    name = "Burned area"
  ) +
  
  scale_fill_manual(
    values = c(
      "rho (MRBA60 vs GFED5)"   = "#75BDD6",
      "rho (MRBA60 vs MCD64A1)" = "#A94444"
    ),
    breaks = c(
      "rho (MRBA60 vs GFED5)",
      "rho (MRBA60 vs MCD64A1)"
    ),
    labels = c(
      "MRBA60 vs GFED5",
      "MRBA60 vs MCD64A1"
    ),
    name = "Correlation"
  ) +
  
  ##--------------------------------------------------------------------------
## Labels
##--------------------------------------------------------------------------
labs(
  x = "Month",
  y = expression("Burned area (km"^2*" month"^{-1}*")")
) +
  
  ##--------------------------------------------------------------------------
## Legend
##--------------------------------------------------------------------------
guides(
  fill = guide_legend(
    order = 1,
    nrow = 1,
    byrow = TRUE,
    override.aes = list(
      alpha = c(0.78, 0.55)
    )
  ),
  color = guide_legend(
    order = 2,
    nrow = 1,
    byrow = TRUE,
    override.aes = list(
      linewidth = 1.05,
      size = 2.4
    )
  ),
  shape = guide_legend(
    order = 2,
    nrow = 1,
    byrow = TRUE
  )
) +
  
  ##--------------------------------------------------------------------------
## Theme
##--------------------------------------------------------------------------
theme_minimal(base_size = 12) +
  
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 13,
      margin = margin(b = 5)
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "center",
    legend.text = element_text(size = 13),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width  = unit(1.05, "cm"),
    legend.spacing.x  = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.20, "cm"),
    legend.margin = margin(t = 8, r = 0, b = 0, l = 0),
    
    panel.grid.major.x = element_line(
      color = "grey86",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_line(
      color = "grey86",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey55",
      fill = NA,
      linewidth = 0.60
    ),
    panel.spacing = unit(0.65, "lines"),
    
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 9)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 9)
    ),
    axis.title.y.right = element_text(
      size = 12,
      margin = margin(l = 9)
    ),
    
    axis.text.x = element_text(
      size = 10,
      color = "grey20"
    ),
    axis.text.y = element_text(
      size = 10,
      color = "grey20"
    ),
    axis.text.y.right = element_text(
      size = 10,
      color = "grey25"
    ),
    
    plot.margin = margin(
      t = 5,
      r = 10,
      b = 14,
      l = 5
    )
  )

print(p)

##------------------------------------------------------------------------------
## Save PDF
##------------------------------------------------------------------------------

ggsave(
  filename = sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.pdf", fig_pdf),
  plot = p,
  device = cairo_pdf,
  width = 10.5,
  height = 13.8,
  units = "in"
)

message(
  "PDF figure saved at: ",
  sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.pdf", fig_pdf)
)

##------------------------------------------------------------------------------
## Save JPEG
##------------------------------------------------------------------------------

ggsave(
  filename = sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.jpeg", fig_pdf),
  plot = p,
  device = "jpeg",
  width = 10.5,
  height = 13.8,
  units = "in",
  dpi = 600,
  bg = "white"
)

message(
  "JPEG figure saved at: ",
  sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.jpeg", fig_pdf)
)

##==============================================================================
## 13B. GLOBAL FIGURE
## - Same concept as regional figure
## - One single global panel
## - Global monthly burned-area climatology
## - Global monthly spatial correlation
##==============================================================================

month_labels <- c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")

##==============================================================================
## 13B.1 GLOBAL MONTHLY BURNED AREA SERIES
##==============================================================================

global_sum_monthly <- function(ba_array, years, months, product_name) {
  
  ba_ts <- sapply(seq_along(years), function(i) {
    sum(ba_array[, , i], na.rm = TRUE)
  })
  
  data.frame(
    product = product_name,
    year    = years,
    month   = months,
    BA      = ba_ts
  )
}

df_global_mrba60 <- global_sum_monthly(
  ba_array     = ba_mrba60[, , 1:ntime, drop = FALSE],
  years        = years,
  months       = months,
  product_name = "MRBA60"
)

df_global_gfed5 <- global_sum_monthly(
  ba_array     = ba_gfed5[, , 1:ntime, drop = FALSE],
  years        = years,
  months       = months,
  product_name = "GFED5"
)

df_global_mcd64a1 <- global_sum_monthly(
  ba_array     = ba_mcd64a1[, , 1:ntime, drop = FALSE],
  years        = years,
  months       = months,
  product_name = "MCD64A1"
)

df_global_all <- dplyr::bind_rows(
  df_global_mrba60,
  df_global_gfed5,
  df_global_mcd64a1
)

clim_global <- df_global_all %>%
  dplyr::group_by(product, month) %>%
  dplyr::summarise(
    BA = mean(BA, na.rm = TRUE),
    .groups = "drop"
  )

##==============================================================================
## 13B.2 GLOBAL SPATIAL SPEARMAN CORRELATION BY MONTH
##==============================================================================

global_monthly_mean_vector <- function(ba_array, idx_m) {
  
  nx <- dim(ba_array)[1]
  ny <- dim(ba_array)[2]
  
  x <- ba_array[, , idx_m, drop = FALSE]
  x <- matrix(x, nrow = nx * ny, ncol = length(idx_m))
  
  rowMeans(x, na.rm = TRUE)
}

compute_global_spatial_corr <- function(ba1, ba2, months_vec, comparison_name, min_n = 30) {
  
  out <- vector("list", 12)
  
  for (m in 1:12) {
    
    idx_m <- which(months_vec == m)
    
    v1 <- global_monthly_mean_vector(ba1, idx_m)
    v2 <- global_monthly_mean_vector(ba2, idx_m)
    
    ok <- is.finite(v1) & is.finite(v2)
    
    n_ok <- sum(ok)
    
    if (n_ok < min_n) {
      rho  <- NA_real_
      pval <- NA_real_
    } else {
      res  <- safe_spearman(v1, v2)
      rho  <- res$rho
      pval <- res$p
      n_ok <- res$n
    }
    
    out[[m]] <- data.frame(
      region_name = "Global",
      month       = m,
      rho         = rho,
      p_value     = pval,
      n_pixels    = n_ok,
      sig_star    = sig_stars(pval),
      significant = is.finite(pval) & pval <= 0.05,
      comparison  = comparison_name
    )
  }
  
  dplyr::bind_rows(out)
}

df_corr_global_gfed <- compute_global_spatial_corr(
  ba1             = ba_mrba60[, , 1:ntime, drop = FALSE],
  ba2             = ba_gfed5[, , 1:ntime, drop = FALSE],
  months_vec      = months,
  comparison_name = "MRBA60 vs GFED5",
  min_n           = 30
)

df_corr_global_mcd <- compute_global_spatial_corr(
  ba1             = ba_mrba60[, , 1:ntime, drop = FALSE],
  ba2             = ba_mcd64a1[, , 1:ntime, drop = FALSE],
  months_vec      = months,
  comparison_name = "MRBA60 vs MCD64A1",
  min_n           = 30
)

df_corr_global_all <- dplyr::bind_rows(
  df_corr_global_gfed,
  df_corr_global_mcd
)

##==============================================================================
## 13B.3 PREPARE GLOBAL BACKGROUND RECTANGLES
##==============================================================================

min_pos_global <- suppressWarnings(
  min(clim_global$BA[is.finite(clim_global$BA) & clim_global$BA > 0], na.rm = TRUE)
)

if (!is.finite(min_pos_global)) {
  stop("There is no BA > 0 in 'clim_global'. Log10 scale cannot be used.")
}

eps_global <- min_pos_global

clim_global_plot <- clim_global %>%
  dplyr::mutate(
    BA_plot = ifelse(is.finite(BA) & BA > 0, BA, eps_global)
  )

rho_min_global <- 0
rho_max_global <- 1.0

logBA_global <- log10(clim_global_plot$BA_plot)

log_min_global <- min(logBA_global[is.finite(logBA_global)], na.rm = TRUE)
log_max_global <- max(logBA_global[is.finite(logBA_global)], na.rm = TRUE)

if (!is.finite(log_min_global) || !is.finite(log_max_global)) {
  stop("Invalid log10(BA_plot) range for global plot.")
}

if (log_max_global == log_min_global) {
  log_min_global <- log_min_global - 1
  log_max_global <- log_max_global + 1
}

a_global <- (log_max_global - log_min_global) / (rho_max_global - rho_min_global)
b_global <- log_min_global - a_global * rho_min_global

df_bg_global_gfed <- df_corr_global_gfed %>%
  dplyr::mutate(
    rho_clip = ifelse(
      is.finite(rho),
      pmin(pmax(rho, rho_min_global), rho_max_global),
      NA_real_
    ),
    ymin = eps_global,
    ymax = ifelse(
      is.finite(rho_clip),
      10^(a_global * rho_clip + b_global),
      NA_real_
    ),
    xmin = month - 0.5,
    xmax = month + 0.5,
    bg_legend = "rho (MRBA60 vs GFED5)",
    ymax = ifelse(is.finite(ymax) & ymax >= ymin, ymax, NA_real_)
  ) %>%
  tibble::as_tibble()

df_bg_global_mcd <- df_corr_global_mcd %>%
  dplyr::mutate(
    rho_clip = ifelse(
      is.finite(rho),
      pmin(pmax(rho, rho_min_global), rho_max_global),
      NA_real_
    ),
    ymin = eps_global,
    ymax = ifelse(
      is.finite(rho_clip),
      10^(a_global * rho_clip + b_global),
      NA_real_
    ),
    xmin = month - 0.5,
    xmax = month + 0.5,
    bg_legend = "rho (MRBA60 vs MCD64A1)",
    ymax = ifelse(is.finite(ymax) & ymax >= ymin, ymax, NA_real_)
  ) %>%
  tibble::as_tibble()

##==============================================================================
## 13B.4 ORDER GLOBAL RECTANGLES
##==============================================================================

df_bg_global_gfed2 <- df_bg_global_gfed %>%
  dplyr::left_join(
    df_bg_global_mcd %>%
      dplyr::select(month, ymax_mcd = ymax),
    by = "month"
  )

df_bg_global_mcd2 <- df_bg_global_mcd %>%
  dplyr::left_join(
    df_bg_global_gfed %>%
      dplyr::select(month, ymax_gfed = ymax),
    by = "month"
  )

df_global_gfed_back <- df_bg_global_gfed2 %>%
  dplyr::filter(
    is.finite(ymax) &
      (is.na(ymax_mcd) | ymax >= ymax_mcd)
  )

df_global_gfed_front <- df_bg_global_gfed2 %>%
  dplyr::filter(
    is.finite(ymax) &
      is.finite(ymax_mcd) &
      ymax < ymax_mcd
  )

df_global_mcd_back <- df_bg_global_mcd2 %>%
  dplyr::filter(
    is.finite(ymax) &
      is.finite(ymax_gfed) &
      ymax > ymax_gfed
  )

df_global_mcd_front <- df_bg_global_mcd2 %>%
  dplyr::filter(
    is.finite(ymax) &
      (is.na(ymax_gfed) | ymax <= ymax_gfed)
  )

##==============================================================================
## 13B.5 GLOBAL PLOT
##==============================================================================

p_global <- ggplot() +
  
  geom_rect(
    data = df_global_gfed_back,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.78,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = df_global_mcd_back,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.55,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = df_global_mcd_front,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.55,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_rect(
    data = df_global_gfed_front,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = bg_legend
    ),
    alpha = 0.78,
    color = NA,
    inherit.aes = FALSE
  ) +
  
  geom_line(
    data = clim_global_plot,
    aes(
      x = month,
      y = BA_plot,
      color = product,
      group = product
    ),
    linewidth = 1.0,
    lineend = "round"
  ) +
  
  geom_point(
    data = clim_global_plot,
    aes(
      x = month,
      y = BA_plot,
      color = product,
      shape = product
    ),
    size = 2.6,
    stroke = 0.25
  ) +
  
  scale_y_log10(
    name = expression("Burned area (km"^2*" month"^{-1}*")"),
    breaks = scales::log_breaks(n = 5),
    labels = scales::label_number(
      big.mark = " ",
      decimal.mark = "."
    ),
    sec.axis = sec_axis(
      transform = ~ (log10(.) - b_global) / a_global,
      breaks = seq(0, 1.0, by = 0.2),
      labels = function(x) sprintf("%.1f", x),
      name = "Spatial correlation"
    )
  ) +
  
  scale_x_continuous(
    breaks = 1:12,
    labels = month_labels,
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_color_manual(
    values = c(
      "MRBA60"  = "#0033FF",
      "GFED5"   = "#202020",
      "MCD64A1" = "#A32020"
    ),
    name = "Burned area"
  ) +
  
  scale_shape_manual(
    values = c(
      "MRBA60"  = 15,
      "GFED5"   = 16,
      "MCD64A1" = 17
    ),
    name = "Burned area"
  ) +
  
  scale_fill_manual(
    values = c(
      "rho (MRBA60 vs GFED5)"   = "#75BDD6",
      "rho (MRBA60 vs MCD64A1)" = "#A94444"
    ),
    breaks = c(
      "rho (MRBA60 vs GFED5)",
      "rho (MRBA60 vs MCD64A1)"
    ),
    labels = c(
      "MRBA60 vs GFED5",
      "MRBA60 vs MCD64A1"
    ),
    name = "Correlation"
  ) +
  
  labs(
    x = "Month",
    y = expression("Burned area (km"^2*" month"^{-1}*")"),
    title = "Global"
  ) +
  
  guides(
    fill = guide_legend(
      order = 1,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        alpha = c(0.78, 0.55)
      )
    ),
    color = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 1.1,
        size = 2.8
      )
    ),
    shape = guide_legend(
      order = 2,
      nrow = 1,
      byrow = TRUE
    )
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16,
      hjust = 0.5,
      margin = margin(b = 8)
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.justification = "center",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width  = unit(1.05, "cm"),
    legend.spacing.x  = unit(0.45, "cm"),
    legend.spacing.y  = unit(0.20, "cm"),
    legend.margin = margin(t = 8, r = 0, b = 0, l = 0),
    
    panel.grid.major.x = element_line(
      color = "grey86",
      linewidth = 0.35
    ),
    panel.grid.major.y = element_line(
      color = "grey86",
      linewidth = 0.35
    ),
    panel.grid.minor = element_blank(),
    
    panel.border = element_rect(
      color = "grey55",
      fill = NA,
      linewidth = 0.60
    ),
    
    axis.title.x = element_text(
      size = 13,
      margin = margin(t = 9)
    ),
    axis.title.y = element_text(
      size = 13,
      margin = margin(r = 9)
    ),
    axis.title.y.right = element_text(
      size = 12,
      margin = margin(l = 9)
    ),
    
    axis.text.x = element_text(
      size = 11,
      color = "grey20"
    ),
    axis.text.y = element_text(
      size = 11,
      color = "grey20"
    ),
    axis.text.y.right = element_text(
      size = 10,
      color = "grey25"
    ),
    
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 14,
      l = 6
    )
  )

print(p_global)

##==============================================================================
## 13B.6 SAVE GLOBAL FIGURE
##==============================================================================

fig_pdf_global <- file.path(out_dir_fig, "Figure9_global.pdf")

ggsave(
  filename = sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.pdf", fig_pdf_global),
  plot = p_global,
  device = cairo_pdf,
  width = 8.5,
  height = 6.5,
  units = "in"
)

ggsave(
  filename = sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.jpeg", fig_pdf_global),
  plot = p_global,
  device = "jpeg",
  width = 8.5,
  height = 6.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

message(
  "Global figure saved at: ",
  sub("\\.pdf$", "_pretty_legend_bottom_no_rho_labels.pdf", fig_pdf_global)
)
##==============================================================================
## 14. SIGNIFICANCE SUMMARY
##==============================================================================

df_sig <- df_corr_all %>%
  dplyr::mutate(
    significant = is.finite(p_value) & p_value <= 0.05,
    month_lab = factor(
      month,
      levels = 1:12,
      labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    )
  )

summary_global <- df_sig %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    n_total  = dplyr::n(),
    n_sig    = sum(significant, na.rm = TRUE),
    n_nonsig = sum(!significant | is.na(significant), na.rm = TRUE),
    all_sig  = n_nonsig == 0,
    .groups = "drop"
  )

by_region <- df_sig %>%
  dplyr::group_by(comparison, region_name) %>%
  dplyr::summarise(
    n_months = dplyr::n(),
    n_sig = sum(significant, na.rm = TRUE),
    n_nonsig = sum(!significant | is.na(significant), na.rm = TRUE),
    pct_sig = 100 * n_sig / n_months,
    .groups = "drop"
  ) %>%
  dplyr::arrange(comparison, n_nonsig, region_name)

print(summary_global)
print(by_region)

##==============================================================================
## 15. PRODUCT-ORDER AND RELATIVE-BIAS SUMMARY
##==============================================================================

clim_w <- clim %>%
  tidyr::pivot_wider(
    names_from = product,
    values_from = BA
  )

req_products <- c("MRBA60", "GFED5", "MCD64A1")

if (!all(req_products %in% names(clim_w))) {
  stop(
    "Some required products are missing in clim_w: ",
    paste(setdiff(req_products, names(clim_w)), collapse = ", ")
  )
}

metrics_order <- clim_w %>%
  dplyr::mutate(
    is_between = MRBA60 >= pmin(GFED5, MCD64A1, na.rm = TRUE) &
      MRBA60 <= pmax(GFED5, MCD64A1, na.rm = TRUE),
    gt_gfed = MRBA60 > GFED5,
    lt_mcd  = MRBA60 < MCD64A1
  )

n_total <- nrow(metrics_order)

summary_order <- metrics_order %>%
  dplyr::summarise(
    n_total = n_total,
    n_between = sum(is_between, na.rm = TRUE),
    pct_between = 100 * n_between / n_total,
    n_gt_gfed = sum(gt_gfed, na.rm = TRUE),
    pct_gt_gfed = 100 * n_gt_gfed / n_total,
    n_lt_mcd = sum(lt_mcd, na.rm = TRUE),
    pct_lt_mcd = 100 * n_lt_mcd / n_total
  )

where_gt_gfed <- metrics_order %>%
  dplyr::filter(gt_gfed) %>%
  dplyr::select(region_name, month, MRBA60, GFED5, MCD64A1) %>%
  dplyr::arrange(region_name, month)

where_lt_mcd <- metrics_order %>%
  dplyr::filter(lt_mcd) %>%
  dplyr::select(region_name, month, MRBA60, GFED5, MCD64A1) %>%
  dplyr::arrange(region_name, month)

bias_rm <- metrics_order %>%
  dplyr::mutate(
    rb_vs_gfed = ifelse(
      is.finite(GFED5) & GFED5 != 0,
      100 * (MRBA60 - GFED5) / GFED5,
      NA_real_
    ),
    rb_vs_mcd = ifelse(
      is.finite(MCD64A1) & MCD64A1 != 0,
      100 * (MRBA60 - MCD64A1) / MCD64A1,
      NA_real_
    )
  )

summary_bias <- bias_rm %>%
  dplyr::summarise(
    mean_rb_vs_gfed   = mean(rb_vs_gfed, na.rm = TRUE),
    median_rb_vs_gfed = median(rb_vs_gfed, na.rm = TRUE),
    mean_rb_vs_mcd    = mean(rb_vs_mcd, na.rm = TRUE),
    median_rb_vs_mcd  = median(rb_vs_mcd, na.rm = TRUE),
    n_rb_gfed = sum(is.finite(rb_vs_gfed)),
    n_rb_mcd  = sum(is.finite(rb_vs_mcd))
  )

summary_rho <- df_corr_all %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    rho_mean = mean(rho, na.rm = TRUE),
    rho_min  = min(rho, na.rm = TRUE),
    rho_max  = max(rho, na.rm = TRUE),
    n_total  = dplyr::n(),
    n_sig    = sum(is.finite(p_value) & p_value <= 0.05, na.rm = TRUE),
    .groups = "drop"
  )

mrba60_clim <- clim %>%
  dplyr::filter(product == "MRBA60") %>%
  dplyr::select(region_name, month, BA_mrba60 = BA) %>%
  dplyr::mutate(
    log10_BA_mrba60 = log10(ifelse(BA_mrba60 > 0, BA_mrba60, NA_real_))
  )

corr_join <- df_corr_all %>%
  dplyr::left_join(
    mrba60_clim,
    by = c("region_name", "month")
  )

spearman_logBA_rho <- corr_join %>%
  dplyr::group_by(comparison) %>%
  dplyr::summarise(
    spearman = suppressWarnings(
      cor(
        log10_BA_mrba60,
        rho,
        method = "spearman",
        use = "complete.obs"
      )
    ),
    n_pairs = sum(is.finite(log10_BA_mrba60) & is.finite(rho)),
    .groups = "drop"
  )

print(summary_order)
print(summary_bias)
print(summary_rho)
print(spearman_logBA_rho)

##==============================================================================
## 16. SAVE TABLES
##==============================================================================

tables_to_save <- list(
  df_mrba60          = df_mrba60,
  df_gfed5           = df_gfed5,
  df_mcd64a1         = df_mcd64a1,
  df_all             = df_all,
  clim               = clim,
  df_corr_gfed       = df_corr_gfed,
  df_corr_mcd        = df_corr_mcd,
  df_corr_all        = df_corr_all,
  df_sig             = df_sig,
  summary_global     = summary_global,
  by_region          = by_region,
  df_bg_gfed         = df_bg_gfed,
  df_bg_mcd          = df_bg_mcd,
  df_gfed_back       = df_gfed_back,
  df_gfed_front      = df_gfed_front,
  df_mcd_back        = df_mcd_back,
  df_mcd_front       = df_mcd_front,
  metrics_order      = metrics_order,
  summary_order      = summary_order,
  where_gt_gfed      = where_gt_gfed,
  where_lt_mcd       = where_lt_mcd,
  bias_rm            = bias_rm,
  summary_bias       = summary_bias,
  summary_rho        = summary_rho,
  corr_join          = corr_join,
  spearman_logBA_rho = spearman_logBA_rho
)

save_tables(
  tables_to_save = tables_to_save,
  out_dir = out_dir_tab
)

message("Tables saved at: ", out_dir_tab)

##==============================================================================
## END
##==============================================================================

gc()
