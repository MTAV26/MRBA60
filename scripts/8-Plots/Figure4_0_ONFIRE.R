# ==========================================================
# LIMPIEZA INICIAL
# ==========================================================
rm(list = ls())
gc()

# ───────────────────────────── LIBRERÍAS ─────────────────────────────
library(ncdf4)
library(lubridate)
library(fields)
library(terra)
library(maps)
library(ggplot2)
library(gridExtra)
library(Metrics)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(sf)
library(rnaturalearth)
library(scales)

# ───────────────────────────── PARÁMETROS ────────────────────────────
nc_dir <- "/mnt/disco6tb/ONFIRE"
files  <- list.files(nc_dir, pattern = "\\.nc$", full.names = TRUE)

# Periodo máximo global de lectura ONFIRE.
# Después se recorta por región:
# Europe      -> 2003–2015
# Canada-NBAC -> 2003–2020
# Resto       -> 2003–2021
start_dt <- as.Date("2003-01-01")
end_dt   <- as.Date("2021-12-31")

# ────────────────────── FUNCIÓN DE LECTURA .NC ONFIRE ───────────────
read_nc_array2 <- function(fpath) {
  
  nc <- nc_open(fpath)
  on.exit(nc_close(nc), add = TRUE)
  
  time_dim <- names(nc$dim)[grepl("time", names(nc$dim), ignore.case = TRUE)][1]
  time_raw <- ncvar_get(nc, time_dim)
  
  units_attr    <- ncatt_get(nc, time_dim, "units")$value
  calendar_attr <- ncatt_get(nc, time_dim, "calendar")$value
  
  if (is.null(calendar_attr) || calendar_attr == "") {
    calendar_attr <- "standard"
  }
  
  if (grepl("^days since", units_attr)) {
    
    origin <- as.Date(sub("days since ", "", units_attr))
    dates  <- origin + time_raw
    
  } else if (grepl("^hours since", units_attr)) {
    
    origin <- as.POSIXct(sub("hours since ", "", units_attr), tz = "UTC")
    dates  <- as.Date(origin + time_raw * 3600)
    
  } else if (grepl("^months since", units_attr)) {
    
    origin <- as.Date(sub("months since ", "", units_attr))
    dates  <- seq(origin, by = "1 month", length.out = length(time_raw))
    
  } else {
    
    stop("Unidades de tiempo no reconocidas: ", units_attr)
  }
  
  sel <- which(dates >= start_dt & dates <= end_dt)
  
  if (length(sel) == 0) {
    stop("No hay datos en ese rango en ", basename(fpath))
  }
  
  varname   <- names(nc$var)[1]
  start_idx <- c(1, 1, sel[1])
  count_idx <- c(-1, -1, length(sel))
  
  arr <- ncvar_get(nc, varname, start = start_idx, count = count_idx)
  
  lon <- ncvar_get(nc, "lon")
  lat <- ncvar_get(nc, "lat")
  
  list(
    file  = basename(fpath),
    data  = arr,
    dates = dates[sel],
    lon   = lon,
    lat   = lat
  )
}

# ─────────────────────── LEER TODOS LOS ARCHIVOS ONFIRE ───────────────
nc_arrays <- lapply(files, read_nc_array2)
names(nc_arrays) <- basename(files)

# ─────────────────────── LEER SHAPEFILE EUROPA ────────────────────────
domain <- vect("/mnt/disco6tb/ONFIRE/domain/domain.shp")

# ─────────────────────── CONVERTIR ONFIRE A km² ───────────────────────
nc_arrays <- lapply(nc_arrays, function(reg) {
  reg$data <- reg$data / 1e6
  reg
})

# ─────────────────────── RECORTE TEMPORAL ESPECÍFICO ──────────────────
# EUROPE:      2003–2015
# CANADA_NBAC: 2003–2020
# Resto:       2003–2021
for (reg in names(nc_arrays)) {
  
  dates_reg <- nc_arrays[[reg]]$dates
  
  if (grepl("EUROPE", reg, ignore.case = TRUE)) {
    
    sel <- which(
      dates_reg >= as.Date("2003-01-01") &
        dates_reg <= as.Date("2015-12-31")
    )
    
  } else if (grepl("CANADA_NBAC", reg, ignore.case = TRUE)) {
    
    sel <- which(
      dates_reg >= as.Date("2003-01-01") &
        dates_reg <= as.Date("2020-12-31")
    )
    
  } else {
    
    sel <- which(
      dates_reg >= as.Date("2003-01-01") &
        dates_reg <= as.Date("2021-12-31")
    )
  }
  
  if (length(sel) == 0) {
    warning("No hay fechas seleccionadas para: ", reg)
    next
  }
  
  nc_arrays[[reg]]$data  <- nc_arrays[[reg]]$data[, , sel, drop = FALSE]
  nc_arrays[[reg]]$dates <- dates_reg[sel]
}

# ─────────────────────── CARGAR FireCCI60 / FireCCI51 / MCD64 / GFED ──
load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/latitude.RData")
load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/longitude.RData")

# Asegurar nombres estándar lon/lat
if (!exists("lon") && exists("longitude")) lon <- longitude
if (!exists("lat") && exists("latitude"))  lat <- latitude

# ─────────────────────── MRBA60 / FireCCI60 ───────────────────────────
Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir       <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv   <- file.path(output_dir, "csv")
output_dir_plot  <- file.path(output_dir, "plot")
output_dir_RData <- file.path(output_dir, "RData")

ruta_RData_MRBA60 <- file.path(
  output_dir_RData,
  "MRBA60_BA_m2_monthly_2003_2024.RData"
)

load(ruta_RData_MRBA60)

BA_final <- BA_MRBA60 / 1e6
BA_final[is.na(BA_final)] <- 0

# ─────────────────────── FireCCI51 2003–2022 ──────────────────────────
load("/mnt/disco6tb/FireCCI60/data_025/out-verifications-2019-2022_025/FireCCI51_2001_2022_0.25degree-download.RData")

BA_FireCCI51 <- Fire51 / 1e6
BA_FireCCI51 <- BA_FireCCI51[, , 25:264]
BA_FireCCI51[is.na(BA_FireCCI51)] <- 0

# ─────────────────────── MCD64A1 2003–2022 ────────────────────────────
nc_path_mcd <- "/mnt/disco6tb/MCD64A1_CMG/nc_output/MCD64CMQ_BA_2001_2022.nc"

nc_mcd   <- nc_open(nc_path_mcd)
BA_MCD64 <- ncvar_get(nc_mcd, "Band1") / 1e4
BA_MCD64 <- BA_MCD64[, , 25:264]
nc_close(nc_mcd)

BA_MCD64[is.na(BA_MCD64)] <- 0

# ─────────────────────── GFED5 2003–2024 ──────────────────────────────
nc_gfed <- nc_open("/mnt/disco6tb/GFED51/burned_area_only/burned_area_2003_2024.nc")

BA_GFED5 <- ncvar_get(nc_gfed, "burned_area") / 1e6

nc_close(nc_gfed)

BA_GFED5[is.na(BA_GFED5)] <- 0

# ─────────────────────── DEFINIR GRILLA 1° ────────────────────────────
lon1 <- seq(
  floor(min(lon)) + 0.5,
  ceiling(max(lon)) - 0.5,
  by = 1
)

lat1 <- seq(
  floor(min(lat)) + 0.5,
  ceiling(max(lat)) - 0.5,
  by = 1
)

# ─────────────────────── FUNCIÓN RE-ASIGNACIÓN 0.25° → 1° ─────────────
regrid_to_1deg <- function(BA_025, lon, lat, lon1, lat1) {
  
  BA_1deg <- array(
    NA_real_,
    dim = c(length(lon1), length(lat1), dim(BA_025)[3])
  )
  
  for (k in seq_len(dim(BA_025)[3])) {
    
    message("Regridding month ", k, " / ", dim(BA_025)[3])
    
    for (i in seq_along(lon1)) {
      
      lon_min <- lon1[i] - 0.5
      lon_max <- lon1[i] + 0.5
      
      idx_lon <- which(lon >= lon_min & lon < lon_max)
      
      if (length(idx_lon) == 0) next
      
      for (j in seq_along(lat1)) {
        
        lat_min <- lat1[j] - 0.5
        lat_max <- lat1[j] + 0.5
        
        idx_lat <- which(lat >= lat_min & lat < lat_max)
        
        if (length(idx_lat) == 0) next
        
        BA_1deg[i, j, k] <- sum(BA_025[idx_lon, idx_lat, k], na.rm = TRUE)
      }
    }
  }
  
  BA_1deg
}

# ─────────────────────── RE-ASIGNAR LOS 4 PRODUCTOS A 1° ──────────────
BA_1deg_FC60  <- regrid_to_1deg(BA_final,     lon, lat, lon1, lat1)
BA_1deg_FC51  <- regrid_to_1deg(BA_FireCCI51, lon, lat, lon1, lat1)
BA_1deg_MCD64 <- regrid_to_1deg(BA_MCD64,     lon, lat, lon1, lat1)
BA_1deg_GFED5 <- regrid_to_1deg(BA_GFED5,     lon, lat, lon1, lat1)

# ─────────────────────── FUNCIÓN RECORTE POR REGIÓN 1° ────────────────
recortar_array_1deg <- function(arr_1deg, lon1, lat1, lon_reg, lat_reg) {
  
  lon_min <- min(lon_reg)
  lon_max <- max(lon_reg)
  lat_min <- min(lat_reg)
  lat_max <- max(lat_reg)
  
  idx_lon <- which(lon1 >= lon_min & lon1 <= lon_max)
  idx_lat <- which(lat1 >= lat_min & lat1 <= lat_max)
  
  arr_recortado <- arr_1deg[idx_lon, idx_lat, , drop = FALSE]
  
  list(
    data = arr_recortado,
    lon  = lon1[idx_lon],
    lat  = lat1[idx_lat]
  )
}

# ─────────────────────── RECORTES POR REGIÓN PARA LOS 4 PRODUCTOS ─────
recortes_FireCCI60 <- list()
recortes_FireCCI51 <- list()
recortes_MCD64     <- list()
recortes_GFED5     <- list()

for (i in seq_along(nc_arrays)) {
  
  nombre_region <- names(nc_arrays)[i]
  
  lon_reg <- nc_arrays[[i]]$lon
  lat_reg <- nc_arrays[[i]]$lat
  
  recortes_FireCCI60[[nombre_region]] <- recortar_array_1deg(
    BA_1deg_FC60, lon1, lat1, lon_reg, lat_reg
  )
  
  recortes_FireCCI51[[nombre_region]] <- recortar_array_1deg(
    BA_1deg_FC51, lon1, lat1, lon_reg, lat_reg
  )
  
  recortes_MCD64[[nombre_region]] <- recortar_array_1deg(
    BA_1deg_MCD64, lon1, lat1, lon_reg, lat_reg
  )
  
  recortes_GFED5[[nombre_region]] <- recortar_array_1deg(
    BA_1deg_GFED5, lon1, lat1, lon_reg, lat_reg
  )
}

# ─────────────────────── FUNCIÓN MÁSCARA VECTORIAL ────────────────────
aplicar_mascara <- function(arr_in, lon_vec, lat_vec, domain) {
  
  if (is.null(arr_in) || is.null(lon_vec) || is.null(lat_vec)) {
    stop("Alguno de los argumentos es NULL")
  }
  
  lon_vec <- as.numeric(lon_vec)
  lat_vec <- as.numeric(lat_vec)
  
  grid_df <- expand.grid(
    lon = lon_vec,
    lat = lat_vec,
    KEEP.OUT.ATTRS = FALSE
  )
  
  coords_mat <- as.matrix(grid_df[, 1:2])
  
  puntos <- terra::vect(
    coords_mat,
    type = "points",
    crs = "EPSG:4326"
  )
  
  dentro <- rowSums(terra::relate(puntos, domain, "intersects")) > 0
  
  mask2d <- matrix(
    dentro,
    nrow = length(lon_vec),
    ncol = length(lat_vec),
    byrow = FALSE
  )
  
  mask3d <- array(mask2d, dim = dim(arr_in))
  
  arr_in[!mask3d] <- NA
  
  arr_in
}

# ─────────────────────── MÁSCARA PARA EUROPA ──────────────────────────
for (reg in names(recortes_FireCCI60)) {
  
  if (!grepl("EUROPE", reg, ignore.case = TRUE)) next
  
  recortes_FireCCI60[[reg]]$data <- aplicar_mascara(
    recortes_FireCCI60[[reg]]$data,
    recortes_FireCCI60[[reg]]$lon,
    recortes_FireCCI60[[reg]]$lat,
    domain
  )
  
  recortes_FireCCI51[[reg]]$data <- aplicar_mascara(
    recortes_FireCCI51[[reg]]$data,
    recortes_FireCCI51[[reg]]$lon,
    recortes_FireCCI51[[reg]]$lat,
    domain
  )
  
  recortes_MCD64[[reg]]$data <- aplicar_mascara(
    recortes_MCD64[[reg]]$data,
    recortes_MCD64[[reg]]$lon,
    recortes_MCD64[[reg]]$lat,
    domain
  )
  
  recortes_GFED5[[reg]]$data <- aplicar_mascara(
    recortes_GFED5[[reg]]$data,
    recortes_GFED5[[reg]]$lon,
    recortes_GFED5[[reg]]$lat,
    domain
  )
}

# ─────────────────────── MÁSCARA PARA OTRAS REGIONES ──────────────────
world_sf <- ne_countries(scale = "medium", returnclass = "sf") %>%
  st_transform(4326)

region_pais_map <- list(
  "AUSTRALIA"   = "Australia",
  "CANADA_NFDB" = "Canada",
  "CANADA_NBAC" = "Canada",
  "CHILE"       = "Chile",
  "US_MTBS"     = "United States of America",
  "US_FPA_FOD"  = "United States of America"
)

for (reg in names(recortes_FireCCI60)) {
  
  if (grepl("EUROPE", reg, ignore.case = TRUE)) next
  
  pais_clave <- names(region_pais_map)[
    sapply(
      names(region_pais_map),
      function(x) grepl(x, reg, ignore.case = TRUE)
    )
  ]
  
  if (length(pais_clave) == 0) {
    warning("No se encontró país para ", reg)
    next
  }
  
  pais_nombre <- region_pais_map[[pais_clave[1]]]
  
  dominio_pais <- world_sf %>%
    dplyr::filter(admin == pais_nombre) %>%
    vect()
  
  recortes_FireCCI60[[reg]]$data <- aplicar_mascara(
    recortes_FireCCI60[[reg]]$data,
    recortes_FireCCI60[[reg]]$lon,
    recortes_FireCCI60[[reg]]$lat,
    dominio_pais
  )
  
  recortes_FireCCI51[[reg]]$data <- aplicar_mascara(
    recortes_FireCCI51[[reg]]$data,
    recortes_FireCCI51[[reg]]$lon,
    recortes_FireCCI51[[reg]]$lat,
    dominio_pais
  )
  
  recortes_MCD64[[reg]]$data <- aplicar_mascara(
    recortes_MCD64[[reg]]$data,
    recortes_MCD64[[reg]]$lon,
    recortes_MCD64[[reg]]$lat,
    dominio_pais
  )
  
  recortes_GFED5[[reg]]$data <- aplicar_mascara(
    recortes_GFED5[[reg]]$data,
    recortes_GFED5[[reg]]$lon,
    recortes_GFED5[[reg]]$lat,
    dominio_pais
  )
}

# ─────────────────────── AJUSTE TEMPORAL PRODUCTOS vs ONFIRE ──────────
for (reg in names(nc_arrays)) {
  
  ntime_loc <- dim(nc_arrays[[reg]]$data)[3]
  
  if (!(reg %in% names(recortes_FireCCI60))) next
  
  dim_p60 <- dim(recortes_FireCCI60[[reg]]$data)[3]
  dim_p51 <- dim(recortes_FireCCI51[[reg]]$data)[3]
  dim_pm  <- dim(recortes_MCD64[[reg]]$data)[3]
  dim_pg  <- dim(recortes_GFED5[[reg]]$data)[3]
  
  if (any(c(dim_p60, dim_p51, dim_pm, dim_pg) < ntime_loc)) {
    
    warning(
      "Menos meses en algún producto que en ONFIRE para región: ", reg,
      " | ONFIRE=", ntime_loc,
      " FC60=", dim_p60,
      " FC51=", dim_p51,
      " MCD64=", dim_pm,
      " GFED=", dim_pg
    )
    
    nmin <- min(c(ntime_loc, dim_p60, dim_p51, dim_pm, dim_pg))
    
    nc_arrays[[reg]]$data  <- nc_arrays[[reg]]$data[, , 1:nmin, drop = FALSE]
    nc_arrays[[reg]]$dates <- nc_arrays[[reg]]$dates[1:nmin]
    
    ntime_loc <- nmin
  }
  
  recortes_FireCCI60[[reg]]$data <- recortes_FireCCI60[[reg]]$data[, , 1:ntime_loc, drop = FALSE]
  recortes_FireCCI51[[reg]]$data <- recortes_FireCCI51[[reg]]$data[, , 1:ntime_loc, drop = FALSE]
  recortes_MCD64[[reg]]$data     <- recortes_MCD64[[reg]]$data[, , 1:ntime_loc, drop = FALSE]
  recortes_GFED5[[reg]]$data     <- recortes_GFED5[[reg]]$data[, , 1:ntime_loc, drop = FALSE]
}

# ─────────────────────── VALIDACIÓN: SCATTER POR REGIÓN ───────────────
dir_plot <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure4_R1/"

if (!dir.exists(dir_plot)) {
  dir.create(dir_plot, recursive = TRUE)
}

col_prod <- c(
  "FireCCI60" = "blue",
  "FireCCI51" = "#E69F00",
  "MCD64A1"   = "brown4",
  "GFED5"     = "gray30"
)

plots_regiones   <- list()
metrics_regiones <- list()

for (reg in names(nc_arrays)) {
  
  message("Procesando región: ", reg)
  
  if (!(reg %in% names(recortes_FireCCI60))) {
    warning("No hay recortes globales para la región: ", reg)
    next
  }
  
  arr_local <- nc_arrays[[reg]]$data
  arr_fc60  <- recortes_FireCCI60[[reg]]$data
  arr_fc51  <- recortes_FireCCI51[[reg]]$data
  arr_mcd   <- recortes_MCD64[[reg]]$data
  arr_gf    <- recortes_GFED5[[reg]]$data
  
  # Comprobación de dimensiones
  if (
    !all(dim(arr_local) == dim(arr_fc60)) ||
    !all(dim(arr_local) == dim(arr_fc51)) ||
    !all(dim(arr_local) == dim(arr_mcd))  ||
    !all(dim(arr_local) == dim(arr_gf))
  ) {
    
    warning("Dimensiones no coinciden en región: ", reg)
    
    print(
      rbind(
        local = dim(arr_local),
        fc60  = dim(arr_fc60),
        fc51  = dim(arr_fc51),
        mcd64 = dim(arr_mcd),
        gfed5 = dim(arr_gf)
      )
    )
    
    next
  }
  
  fechas_reg <- nc_arrays[[reg]]$dates
  anios_reg  <- format(fechas_reg, "%Y")
  anios_u    <- sort(unique(anios_reg))
  nA         <- length(anios_u)
  
  sum_local <- numeric(nA)
  sum_fc60  <- numeric(nA)
  sum_fc51  <- numeric(nA)
  sum_mcd64 <- numeric(nA)
  sum_gfed5 <- numeric(nA)
  
  for (k in seq_len(nA)) {
    
    idx_t <- which(anios_reg == anios_u[k])
    
    sum_local[k] <- sum(arr_local[, , idx_t], na.rm = TRUE)
    sum_fc60[k]  <- sum(arr_fc60[, , idx_t],  na.rm = TRUE)
    sum_fc51[k]  <- sum(arr_fc51[, , idx_t],  na.rm = TRUE)
    sum_mcd64[k] <- sum(arr_mcd[, , idx_t],   na.rm = TRUE)
    sum_gfed5[k] <- sum(arr_gf[, , idx_t],    na.rm = TRUE)
  }
  
  df_reg <- data.frame(
    year  = as.integer(anios_u),
    local = sum_local,
    fc60  = sum_fc60,
    fc51  = sum_fc51,
    mcd64 = sum_mcd64,
    gfed5 = sum_gfed5
  )
  
  df_reg_long <- df_reg %>%
    pivot_longer(
      cols      = c(fc60, fc51, mcd64, gfed5),
      names_to  = "dataset",
      values_to = "value"
    ) %>%
    mutate(
      dataset = factor(
        dataset,
        levels = c("fc60", "fc51", "mcd64", "gfed5"),
        labels = c("FireCCI60", "FireCCI51", "MCD64A1", "GFED5")
      ),
      region = reg
    )
  
  # ─────────────── MÉTRICAS ───────────────
  grupos_reg <- df_reg_long %>%
    group_split(dataset, .keep = TRUE)
  
  metrics_reg <- purrr::map_dfr(grupos_reg, function(d) {
    
    ds <- as.character(unique(d$dataset))
    
    if (nrow(d) == 0) {
      
      return(
        tibble(
          region       = reg,
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
        )
      )
    }
    
    fit <- lm(value ~ local, data = d)
    
    a <- unname(coef(fit)[1])
    b <- unname(coef(fit)[2])
    
    p_slope <- summary(fit)$coefficients[2, 4]
    
    err <- d$value - d$local
    
    rel_bias <- err / d$local
    rel_bias[d$local == 0] <- NA_real_
    
    bias_pct_val <- mean(rel_bias, na.rm = TRUE) * 100
    
    sp <- tryCatch({
      
      ok <- is.finite(d$local) & is.finite(d$value)
      
      x <- d$local[ok]
      y <- d$value[ok]
      
      if (length(x) < 3 || sd(x) == 0 || sd(y) == 0) {
        
        list(rho = NA_real_, p = NA_real_)
        
      } else {
        
        ct <- suppressWarnings(
          cor.test(x, y, method = "spearman", exact = FALSE)
        )
        
        list(
          rho = unname(ct$estimate),
          p   = unname(ct$p.value)
        )
      }
      
    }, error = function(e) {
      
      list(rho = NA_real_, p = NA_real_)
    })
    
    tibble(
      region       = reg,
      dataset      = ds,
      r2           = summary(fit)$r.squared,
      r2sig        = ifelse(!is.na(p_slope) && p_slope < 0.05, "*", ""),
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
    mutate(
      dataset = factor(
        dataset,
        levels = c("FireCCI60", "FireCCI51", "MCD64A1", "GFED5")
      )
    ) %>%
    arrange(dataset)
  
  metrics_regiones[[reg]] <- metrics_reg
  
  # ─────────────── TABLA CSV POR REGIÓN ───────────────
  tabla_stats_reg <- metrics_reg %>%
    transmute(
      Región         = region,
      Producto       = as.character(dataset),
      `Ecuación`     = paste0(
        "y = ",
        round(intercept, 3),
        " + ",
        round(slope, 3),
        "·x"
      ),
      `R²`           = paste0(round(r2, 3), r2sig),
      `RMSE (km²)`   = round(rmse, 3),
      `MAE (km²)`    = round(mae, 3),
      `BIAS (km²)`   = round(bias, 3),
      `BIAS (%)`     = round(bias_pct, 2),
      `Spearman (ρ)` = paste0(round(spearman, 3), spearman_sig),
      `p Spearman`   = signif(spearman_p, 3)
    )
  
  nombre_base_csv <- gsub("\\.nc$", "", reg)
  
  write.csv(
    tabla_stats_reg,
    file = file.path(
      dir_plot,
      paste0("validacion_ONFIRE_", nombre_base_csv, "_stats.csv")
    ),
    row.names = FALSE
  )
  
  # ─────────────── TÍTULO Y LETRAS ───────────────
  nombre_base <- gsub("^Monthly_Burned_area_\\d{4}_\\d{4}_", "", reg)
  nombre_base <- gsub("_v1\\.nc$", "", nombre_base)
  
  nombre_base <- gsub("CHILE", "Chile", nombre_base)
  nombre_base <- gsub("EUROPE", "Europe", nombre_base)
  nombre_base <- gsub("AUSTRALIA", "Australia", nombre_base)
  nombre_base <- gsub("CANADA_NFDB", "Canada_NFDB", nombre_base)
  nombre_base <- gsub("CANADA_NBAC", "Canada-NBAC", nombre_base)
  nombre_base <- gsub("US_FPA_FOD", "US_FPA_FOD", nombre_base)
  nombre_base <- gsub("US_MTBS", "US_MTBS", nombre_base)
  
  letras_map <- list(
    "Chile"       = "a",
    "Canada_NFDB" = "b",
    "US_FPA_FOD"  = "c",
    "US_MTBS"     = "d",
    "Australia"   = "e",
    "Canada-NBAC" = "f",
    "Europe"      = "g"
  )
  
  periodo_texto <- dplyr::case_when(
    nombre_base == "Europe"      ~ "2003–2015",
    nombre_base == "Canada-NBAC" ~ "2003–2020",
    TRUE                         ~ "2003–2021"
  )
  
  letra <- letras_map[[nombre_base]]
  
  if (is.null(letra)) {
    letra <- ""
  }
  
  titulo_reg <- paste0(letra, ") ONFIRE (", nombre_base, "; ", periodo_texto, ")")
  
  print(titulo_reg)
  
  # ─────────────── SCATTER ───────────────
  p_reg <- ggplot(
    df_reg_long,
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
    scale_color_manual(values = col_prod) +
    scale_fill_manual(values = col_prod) +
    scale_x_continuous(
      labels = scales::number_format(accuracy = 1)
    ) +
    scale_y_continuous(
      labels = scales::number_format(accuracy = 1)
    ) +
    labs(
      title = titulo_reg,
      x = expression("ONFIRE BA (km"^2~yr^{-1}*")"),
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
  
  print(p_reg)
  
  plots_regiones[[reg]] <- p_reg
  
  ggsave(
    filename = file.path(
      dir_plot,
      paste0(
        "validacion_ONFIRE_",
        nombre_base,
        "_",
        min(df_reg$year),
        "_",
        max(df_reg$year),
        "_4prod.pdf"
      )
    ),
    plot   = p_reg,
    width  = 12,
    height = 7,
    dpi    = 300,
    device = "pdf"
  )
  
  ggsave(
    filename = file.path(
      dir_plot,
      paste0(
        "validacion_ONFIRE_",
        nombre_base,
        "_",
        min(df_reg$year),
        "_",
        max(df_reg$year),
        "_4prod.jpeg"
      )
    ),
    plot   = p_reg,
    width  = 12,
    height = 7,
    dpi    = 300,
    device = "jpeg"
  )
}

# ─────────────────────── TABLA GLOBAL TODAS LAS REGIONES ──────────────
tabla_stats_onfire <- bind_rows(metrics_regiones) %>%
  mutate(
    r2           = round(r2, 3),
    slope        = round(slope, 3),
    intercept    = round(intercept, 3),
    rmse         = round(rmse, 3),
    mae          = round(mae, 3),
    bias         = round(bias, 3),
    bias_pct     = round(bias_pct, 2),
    spearman     = round(spearman, 3),
    spearman_p   = signif(spearman_p, 3),
    spearman_sig = ifelse(!is.na(spearman_p) & spearman_p < 0.05, "*", "")
  )

write.csv(
  tabla_stats_onfire,
  file = file.path(
    dir_plot,
    "validacion_ONFIRE_todas_regiones_stats_4prod_periodos_region.csv"
  ),
  row.names = FALSE
)

# ─────────────────────── MOSTRAR EUROPE SI EXISTE ─────────────────────
idx_eu_plot <- grep("EUROPE", names(plots_regiones), ignore.case = TRUE)

if (length(idx_eu_plot) > 0) {
  print(plots_regiones[[idx_eu_plot[1]]])
}