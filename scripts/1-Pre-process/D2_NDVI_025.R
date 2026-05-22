rm(list = ls())
graphics.off()
gc()

library(terra)
library(ncdf4)

# ============================================================
# 1) Rutas
# ============================================================
dir_in   <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/MOD13C2_NDVI"
dir_out  <- "/mnt/disco6tb/MRBA60/data/A2_TEMP"
out_nc   <- file.path(dir_out, "MOD13C2_NDVI_2003_2024_025.nc")

file_lon <- file.path(dir_out, "longitude.RData")
file_lat <- file.path(dir_out, "latitude.RData")

dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2) Parámetros
# ============================================================
apply_scale  <- TRUE
scale_factor <- 0.0001
fillvalue    <- 1e20

# TRUE = media ponderada por área (recomendado)
# FALSE = media simple 5x5
use_area_weighted_mean <- TRUE

# tolerancia para comparar geometrías
tol_res <- 1e-10
tol_ext <- 1e-8

# ============================================================
# 3) Funciones auxiliares
# ============================================================
get_date_from_filename <- function(fname) {
  posA <- regexpr("A\\d{7}", fname)
  if (posA[1] == -1) return(as.Date(NA))
  
  adate <- regmatches(fname, posA)   # "A2000032"
  year  <- as.integer(substr(adate, 2, 5))
  doy   <- as.integer(substr(adate, 6, 8))
  
  as.Date(doy - 1, origin = paste0(year, "-01-01"))
}

load_numeric_vector_from_rdata <- function(path) {
  e <- new.env()
  objs <- load(path, envir = e)
  if (length(objs) == 0) stop("No hay objetos en: ", path)
  
  # prioriza vectores numéricos
  idx <- sapply(objs, function(nm) {
    x <- e[[nm]]
    is.numeric(x) && is.atomic(x)
  })
  
  if (!any(idx)) {
    stop("No se encontró ningún vector numérico en: ", path)
  }
  
  x <- e[[objs[which(idx)[1]]]]
  as.vector(x)
}

# ============================================================
# 4) Cargar lon/lat objetivo a 0.25°
# ============================================================
lon_out <- load_numeric_vector_from_rdata(file_lon)
lat_out <- load_numeric_vector_from_rdata(file_lat)

lon_out <- as.numeric(lon_out)
lat_out <- as.numeric(lat_out)

# ordenar por seguridad
lon_out <- sort(unique(lon_out))
lat_out <- sort(unique(lat_out))

if (length(lon_out) < 2 || length(lat_out) < 2) {
  stop("Los vectores lon/lat no son válidos.")
}

dx_out <- median(diff(lon_out))
dy_out <- median(diff(lat_out))

if (any(abs(diff(lon_out) - dx_out) > 1e-8)) {
  stop("longitude.RData no parece una malla regular.")
}
if (any(abs(diff(lat_out) - dy_out) > 1e-8)) {
  stop("latitude.RData no parece una malla regular.")
}

# Plantilla objetivo 0.25°
r_template <- rast(
  xmin = min(lon_out) - dx_out/2,
  xmax = max(lon_out) + dx_out/2,
  ymin = min(lat_out) - dy_out/2,
  ymax = max(lat_out) + dy_out/2,
  ncols = length(lon_out),
  nrows = length(lat_out),
  crs  = "EPSG:4326"
)

nlon_out <- ncol(r_template)
nlat_out <- nrow(r_template)

cat("Grid objetivo 0.25°:\n")
cat("nlon =", nlon_out, " | nlat =", nlat_out, "\n")
cat("Lon:", min(lon_out), "->", max(lon_out), "\n")
cat("Lat:", min(lat_out), "->", max(lat_out), "\n")
cat("Res:", dx_out, "x", dy_out, "\n\n")

# ============================================================
# 5) Listar archivos y filtrar 2003-2024
# ============================================================
files <- list.files(
  dir_in,
  pattern = "^MOD13C2\\.A\\d{7}.*\\.nc$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No se encontraron archivos .nc en: ", dir_in)
}

dates <- as.Date(
  sapply(basename(files), get_date_from_filename),
  origin = "1970-01-01"
)

years <- as.integer(format(dates, "%Y"))
sel   <- !is.na(dates) & years >= 2003 & years <= 2024

files <- files[sel]
dates <- dates[sel]

if (length(files) == 0) {
  stop("No hay archivos entre 2003 y 2024.")
}

ord   <- order(dates)
files <- files[ord]
dates <- dates[ord]

if (any(duplicated(dates))) {
  print(data.frame(file = basename(files), date = dates))
  stop("Se detectaron fechas duplicadas.")
}

expected_months <- format(seq(as.Date("2003-01-01"), as.Date("2024-12-01"), by = "month"), "%Y-%m")
found_months    <- format(dates, "%Y-%m")
missing_months  <- setdiff(expected_months, found_months)

if (length(missing_months) > 0) {
  warning("Faltan meses en la serie: ", paste(missing_months, collapse = ", "))
}

cat("Número de archivos seleccionados:", length(files), "\n")
cat("Rango temporal:", as.character(min(dates)), "->", as.character(max(dates)), "\n\n")

# ============================================================
# 6) Raster de referencia (0.05°)
# ============================================================
r0 <- rast(files[1])
if (nlyr(r0) > 1) r0 <- r0[[1]]

nlon_in <- ncol(r0)
nlat_in <- nrow(r0)
res0    <- res(r0)
ext0    <- as.vector(ext(r0))
crs0    <- crs(r0)

cat("Grid original:\n")
cat("nlon =", nlon_in, " | nlat =", nlat_in, "\n")
cat("Res:", res0[1], "x", res0[2], "\n")
cat("CRS:", crs0, "\n\n")

# factor de agregación esperado
fact_x <- dx_out / res0[1]
fact_y <- dy_out / res0[2]

if (abs(fact_x - round(fact_x)) > 1e-8 || abs(fact_y - round(fact_y)) > 1e-8) {
  stop("La relación entre la resolución original y la de salida no es entera.")
}

fact_x <- as.integer(round(fact_x))
fact_y <- as.integer(round(fact_y))

cat("Factor de agregación:", fact_x, "x", fact_y, "\n\n")

# área de celda para media ponderada
if (use_area_weighted_mean) {
  w_area <- cellSize(r0, unit = "km")
}

# ============================================================
# 7) Definir NetCDF de salida
#    Orden final: NDVI(lon, lat, time)
# ============================================================
time_vals <- as.numeric(dates - as.Date("1970-01-01"))

dim_lon <- ncdim_def(
  name = "lon",
  units = "degrees_east",
  vals = lon_out,
  create_dimvar = TRUE
)

dim_lat <- ncdim_def(
  name = "lat",
  units = "degrees_north",
  vals = lat_out,
  create_dimvar = TRUE
)

dim_time <- ncdim_def(
  name = "time",
  units = "days since 1970-01-01 00:00:00",
  vals = time_vals,
  unlim = TRUE,
  create_dimvar = TRUE
)

var_longname <- if (use_area_weighted_mean) {
  "MOD13C2 monthly NDVI aggregated to 0.25 degree using area-weighted mean"
} else {
  "MOD13C2 monthly NDVI aggregated to 0.25 degree using simple mean"
}

var_ndvi <- ncvar_def(
  name = "NDVI",
  units = "1",
  dim = list(dim_lon, dim_lat, dim_time),
  missval = fillvalue,
  longname = var_longname,
  prec = "float",
  compression = 4
)

if (file.exists(out_nc)) file.remove(out_nc)

nc <- nc_create(
  filename = out_nc,
  vars = list(var_ndvi),
  force_v4 = TRUE
)

ncatt_put(nc, "lon",  "standard_name", "longitude")
ncatt_put(nc, "lon",  "axis", "X")
ncatt_put(nc, "lat",  "standard_name", "latitude")
ncatt_put(nc, "lat",  "axis", "Y")
ncatt_put(nc, "time", "standard_name", "time")
ncatt_put(nc, "time", "calendar", "standard")

ncatt_put(nc, "NDVI", "coordinates", "lon lat")
ncatt_put(nc, "NDVI", "scale_applied", ifelse(apply_scale, scale_factor, 1))
ncatt_put(nc, "NDVI", "aggregation_method",
          ifelse(use_area_weighted_mean, "area_weighted_mean", "simple_mean"))

ncatt_put(nc, 0, "title", "MOD13C2 NDVI 2003-2024 aggregated to 0.25 degree")
ncatt_put(nc, 0, "source", "Merged from /mnt/disco6tb/MRBA60/data/A2_TEMP/MOD13C2_NDVI")
ncatt_put(nc, 0, "history", paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                  "- aggregated to 0.25 degree and merged in R"))
ncatt_put(nc, 0, "Conventions", "CF-1.8")
ncatt_put(nc, 0, "crs", ifelse(is.na(crs0) || crs0 == "", "undefined", crs0))

# ============================================================
# 8) Procesar capa a capa y escribir
# ============================================================
for (i in seq_along(files)) {
  cat("Procesando", i, "de", length(files), ":", basename(files[i]), "\n")
  
  r <- rast(files[i])
  if (nlyr(r) > 1) r <- r[[1]]
  
  # ----------------------------
  # Validaciones de grid original
  # ----------------------------
  if (ncol(r) != nlon_in || nrow(r) != nlat_in) {
    nc_close(nc)
    stop("Dimensiones distintas en: ", basename(files[i]))
  }
  
  if (any(abs(res(r) - res0) > tol_res)) {
    nc_close(nc)
    stop("Resolución distinta en: ", basename(files[i]))
  }
  
  if (any(abs(as.vector(ext(r)) - ext0) > tol_ext)) {
    nc_close(nc)
    stop("Extensión distinta en: ", basename(files[i]))
  }
  
  crs_i <- crs(r)
  same_crs <- identical(crs_i, crs0) || (is.na(crs_i) && is.na(crs0))
  if (!same_crs) {
    nc_close(nc)
    stop("CRS distinto en: ", basename(files[i]))
  }
  
  # Escalado
  if (apply_scale) {
    r <- r * scale_factor
  }
  
  # ----------------------------
  # Agregación 0.05 -> 0.25
  # ----------------------------
  if (use_area_weighted_mean) {
    # numerador: sum(NDVI * area)
    num <- aggregate(r * w_area, fact = c(fact_x, fact_y), fun = sum, na.rm = TRUE)
    
    # denominador: sum(area) solo donde NDVI no es NA
    den <- aggregate(mask(w_area, r), fact = c(fact_x, fact_y), fun = sum, na.rm = TRUE)
    
    r025 <- num / den
    r025 <- mask(r025, den)
  } else {
    r025 <- aggregate(r, fact = c(fact_x, fact_y), fun = mean, na.rm = TRUE)
  }
  
  # ----------------------------
  # Comprobar que coincide con la grid objetivo
  # ----------------------------
  if (ncol(r025) != nlon_out || nrow(r025) != nlat_out) {
    nc_close(nc)
    stop("La malla agregada no coincide con la malla objetivo en: ", basename(files[i]))
  }
  
  if (any(abs(res(r025) - c(dx_out, dy_out)) > tol_res)) {
    nc_close(nc)
    stop("La resolución agregada no coincide con 0.25° en: ", basename(files[i]))
  }
  
  ext025 <- as.vector(ext(r025))
  extT   <- as.vector(ext(r_template))
  if (any(abs(ext025 - extT) > tol_ext)) {
    nc_close(nc)
    stop("La extensión agregada no coincide con la malla objetivo en: ", basename(files[i]))
  }
  
  # ----------------------------
  # Extraer valores y escribir como [lon, lat]
  # terra: filas norte->sur, columnas oeste->este
  # netcdf final: lon asc, lat asc
  # ----------------------------
  v <- values(r025, mat = FALSE)
  m <- matrix(v, nrow = nlat_out, ncol = nlon_out, byrow = TRUE)
  m <- m[nlat_out:1, , drop = FALSE]  # sur -> norte
  m <- t(m)                           # [lon, lat]
  
  m[is.na(m)] <- fillvalue
  
  ncvar_put(
    nc,
    varid = var_ndvi,
    vals = m,
    start = c(1, 1, i),
    count = c(nlon_out, nlat_out, 1)
  )
}

nc_close(nc)

cat("\nArchivo final creado:\n", out_nc, "\n")
cat("nlon =", nlon_out, "| nlat =", nlat_out, "| ntime =", length(dates), "\n")
cat("Rango temporal:", as.character(min(dates)), "->", as.character(max(dates)), "\n")
cat("Lon:", min(lon_out), "->", max(lon_out), "\n")
cat("Lat:", min(lat_out), "->", max(lat_out), "\n")
cat("Método de agregación:", ifelse(use_area_weighted_mean,
                                    "media ponderada por área",
                                    "media simple"), "\n")