rm(list = ls())
graphics.off()
gc()

# ==========================================================
# 1) Cargar datos MRBA60
# ==========================================================

load("/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_BA_m2_monthly_2003_2024.RData")
dim(BA_MRBA60)

load("/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/MRBA60_Unc_m2_monthly_2003_2024.RData")
dim(Unc_MRBA60)

burned_area <- BA_MRBA60
uncertainty <- Unc_MRBA60

burned_area[is.na(burned_area)] <- 0
uncertainty[is.na(uncertainty)] <- 0

if (exists("BA_FireCCI60"))  rm(BA_FireCCI60)
if (exists("Unc_FireCCI60")) rm(Unc_FireCCI60)
gc()

# ==========================================================
# 2) Librerías
# ==========================================================

library(ncdf4)

if (!requireNamespace("uuid", quietly = TRUE)) {
  gen_tracking_id <- function() {
    sprintf("tmp-%s", format(Sys.time(), "%Y%m%d%H%M%S"))
  }
} else {
  gen_tracking_id <- function() {
    uuid::UUIDgenerate()
  }
}

# ==========================================================
# 3) Coordenadas y fechas
# ==========================================================

lon <- seq(-180 + 0.125, 180 - 0.125, by = 0.25)  # 1440
lat <- seq(-90  + 0.125,  90 - 0.125, by = 0.25)  # 720

dates <- seq(
  as.Date("2003-01-01"),
  as.Date("2024-12-01"),
  by = "1 month"
)

stopifnot(length(lon) == dim(burned_area)[1])
stopifnot(length(lat) == dim(burned_area)[2])
stopifnot(length(dates) == dim(burned_area)[3])
stopifnot(all(dim(burned_area) == dim(uncertainty)))

half <- 0.125

month_bounds <- function(d) {
  t0 <- d
  t1 <- seq(d, by = "month", length.out = 2)[2] - 1
  c(t0, t1)
}

# ==========================================================
# 4) Carpeta de salida
# ==========================================================

out_dir <- "/mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==========================================================
# 5) Función para escribir un NetCDF mensual
# ==========================================================

write_month_nc <- function(k) {
  
  stopifnot(k >= 1, k <= length(dates))
  
  d <- dates[k]
  year_str <- strftime(d, "%Y")
  
  # ----------------------------------------------------------
  # Nombre de archivo definitivo
  # Ejemplo:
  # 20140101-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc
  # ----------------------------------------------------------
  
  date_str <- strftime(d, "%Y%m01")
  
  out_name <- sprintf(
    "%s-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc",
    date_str
  )
  
  out_dir_year <- file.path(out_dir, year_str)
  
  if (!dir.exists(out_dir_year)) {
    dir.create(out_dir_year, recursive = TRUE)
  }
  
  out_path <- file.path(out_dir_year, out_name)
  
  if (!dir.exists(dirname(out_path))) {
    dir.create(dirname(out_path), recursive = TRUE)
  }
  
  # ==========================================================
  # 5.1) Definir dimensiones
  # ==========================================================
  
  dim_lon <- ncdim_def(
    name = "lon",
    units = "degree_east",
    vals = lon,
    longname = "longitude",
    create_dimvar = TRUE
  )
  
  dim_lat <- ncdim_def(
    name = "lat",
    units = "degree_north",
    vals = lat,
    longname = "latitude",
    create_dimvar = TRUE
  )
  
  dim_time <- ncdim_def(
    name = "time",
    units = "days since 1970-01-01 00:00:00",
    vals = as.numeric(d - as.Date("1970-01-01")),
    longname = "time",
    calendar = "standard",
    unlim = TRUE,
    create_dimvar = TRUE
  )
  
  dim_bounds <- ncdim_def(
    name = "bounds",
    units = "",
    vals = 1:2,
    create_dimvar = FALSE
  )
  
  # ==========================================================
  # 5.2) Definir variables auxiliares
  # ==========================================================
  
  v_lonb <- ncvar_def(
    name = "lon_bounds",
    units = "",
    dim = list(dim_bounds, dim_lon),
    longname = "",
    prec = "double"
  )
  
  v_latb <- ncvar_def(
    name = "lat_bounds",
    units = "",
    dim = list(dim_bounds, dim_lat),
    longname = "",
    prec = "double"
  )
  
  v_timeb <- ncvar_def(
    name = "time_bounds",
    units = "days since 1970-01-01 00:00:00",
    dim = list(dim_bounds, dim_time),
    longname = "",
    prec = "double"
  )
  
  v_crs <- ncvar_def(
    name = "crs",
    units = "",
    dim = list(),
    prec = "integer"
  )
  
  # ==========================================================
  # 5.3) Definir variables de datos
  # ==========================================================
  
  vb <- ncvar_def(
    name = "burned_area",
    units = "m2",
    dim = list(dim_lon, dim_lat, dim_time),
    missval = 9.96921e36,
    longname = "total burned_area",
    prec = "float",
    compression = 5,
    shuffle = FALSE,
    chunksizes = c(length(lon), length(lat), 1)
  )
  
  vu <- ncvar_def(
    name = "uncertainty",
    units = "m2",
    dim = list(dim_lon, dim_lat, dim_time),
    missval = 9.96921e36,
    longname = "uncertainty of burned area, calculated as Root Mean Square Error.",
    prec = "float",
    compression = 5,
    shuffle = FALSE,
    chunksizes = c(length(lon), length(lat), 1)
  )
  
  # ==========================================================
  # 5.4) Crear archivo NetCDF
  # ==========================================================
  
  ncout <- nc_create(
    filename = out_path,
    vars = list(v_lonb, v_latb, v_timeb, v_crs, vb, vu),
    force_v4 = TRUE
  )
  
  # ==========================================================
  # 5.5) Atributos de coordenadas
  # ==========================================================
  
  ncatt_put(ncout, "lon", "standard_name", "longitude")
  ncatt_put(ncout, "lon", "long_name",     "longitude")
  ncatt_put(ncout, "lon", "bounds",        "lon_bounds")
  
  ncatt_put(ncout, "lat", "standard_name", "latitude")
  ncatt_put(ncout, "lat", "long_name",     "latitude")
  ncatt_put(ncout, "lat", "bounds",        "lat_bounds")
  
  ncatt_put(ncout, "time", "standard_name", "time")
  ncatt_put(ncout, "time", "long_name",     "time")
  ncatt_put(ncout, "time", "bounds",        "time_bounds")
  ncatt_put(ncout, "time", "calendar",      "standard")
  
  # ==========================================================
  # 5.6) Atributos CRS
  # ==========================================================
  
  ncatt_put(
    ncout,
    "crs",
    "wkt",
    'GEOGCS["WGS84(DD)",DATUM["WGS84",SPHEROID["WGS84", 6378137.0, 298.257223563]],PRIMEM["Greenwich", 0.0],UNIT["degree", 0.017453292519943295],AXIS["Geodetic longitude", EAST],AXIS["Geodetic latitude", NORTH]]'
  )
  
  ncatt_put(
    ncout,
    "crs",
    "i2m",
    "0.25,0.0,0.0,-0.25,-180.0,90.0"
  )
  
  # ==========================================================
  # 5.7) Atributos de burned_area
  # ==========================================================
  
  ncatt_put(ncout, "burned_area", "standard_name", "burned_area")
  ncatt_put(ncout, "burned_area", "long_name",     "total burned_area")
  ncatt_put(ncout, "burned_area", "units",         "m2")
  ncatt_put(ncout, "burned_area", "cell_methods",  "time: sum")
  ncatt_put(ncout, "burned_area", "valid_range",   as.double(c(0, 769288944)))
  
  # ==========================================================
  # 5.8) Atributos de uncertainty
  # ==========================================================
  
  ncatt_put(ncout, "uncertainty", "units",       "m2")
  ncatt_put(ncout, "uncertainty", "long_name",   "uncertainty of burned area, calculated as Root Mean Square Error.")
  ncatt_put(ncout, "uncertainty", "valid_range", as.double(c(0, 769288944)))
  
  # ==========================================================
  # 5.9) Atributos globales definitivos
  # ==========================================================
  
  tb <- month_bounds(d)
  
  time_coverage_start_val <- strftime(
    tb[1],
    "%Y%m%dT000000Z",
    tz = "UTC"
  )
  
  time_coverage_end_val <- strftime(
    tb[2],
    "%Y%m%dT235959Z",
    tz = "UTC"
  )
  
  date_created_utc <- strftime(
    Sys.time(),
    "%Y%m%dT%H%M%SZ",
    tz = "UTC"
  )
  
  ncatt_put(
    ncout,
    0,
    "title",
    "Harmonised Medium Resolution Burned Area Grid product, version 6.0 (MRBA60)"
  )
  
  ncatt_put(
    ncout,
    0,
    "institution",
    "University of Alcala"
  )
  
  ncatt_put(
    ncout,
    0,
    "source",
    "Harmonised FireCCI51 to FireCCIS311 (now renamed MRBA60) burned area, MODIS MCD14ML thermal anomalies. Ancillary datasets: climatological and vegetation variables (see ATBD of the product for more information)"
  )
  
  ncatt_put(
    ncout,
    0,
    "history",
    paste("Created on", strftime(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  )
  
  ncatt_put(
    ncout,
    0,
    "references",
    "See https://climate.esa.int/en/projects/fire/"
  )
  
  ncatt_put(
    ncout,
    0,
    "tracking_id",
    gen_tracking_id()
  )
  
  ncatt_put(
    ncout,
    0,
    "Conventions",
    "CF-1.7"
  )
  
  ncatt_put(
    ncout,
    0,
    "product_version",
    "v6.0.0"
  )
  
  ncatt_put(
    ncout,
    0,
    "format_version",
    "CCI Data Standards v2.3"
  )
  
  ncatt_put(
    ncout,
    0,
    "summary",
    "MRBA60H is the harmonisation of two burned area datasets: FireCCI51 and MRBA60 (previously known as FireCCIS311), using a model to estimate the burned area detection of MRBA60 for the years when only the FireCCI51 dataset was available (2003-2018)."
  )
  
  ncatt_put(
    ncout,
    0,
    "keywords",
    "Burned Area, Fire Disturbance, Climate Change, ESA, GCOS"
  )
  
  ncatt_put(
    ncout,
    0,
    "id",
    out_name
  )
  
  ncatt_put(
    ncout,
    0,
    "naming_authority",
    "int.esa.climate"
  )
  
  ncatt_put(
    ncout,
    0,
    "doi",
    "10.5285/db75c5f51ee240ae8743355dcebbb9b9"
  )
  
  ncatt_put(
    ncout,
    0,
    "keywords_vocabulary",
    "none"
  )
  
  ncatt_put(
    ncout,
    0,
    "cdm_data_type",
    "Grid"
  )
  
  ncatt_put(
    ncout,
    0,
    "comment",
    "These data were produced as part of the Climate Change Initiative Programme, Fire Disturbance ECV."
  )
  
  ncatt_put(
    ncout,
    0,
    "date_created",
    date_created_utc
  )
  
  ncatt_put(
    ncout,
    0,
    "creator_name",
    "University of Alcala"
  )
  
  ncatt_put(
    ncout,
    0,
    "creator_url",
    "https://geogra.uah.es/gita/en/"
  )
  
  ncatt_put(
    ncout,
    0,
    "creator_email",
    "emilio.chuvieco@uah.es"
  )
  
  ncatt_put(
    ncout,
    0,
    "contact",
    "mlucrecia.pettinari@uah.es"
  )
  
  ncatt_put(
    ncout,
    0,
    "developer_email",
    "miguela.torres@uah.es"
  )
  
  ncatt_put(
    ncout,
    0,
    "project",
    "Climate Change Initiative - European Space Agency"
  )
  
  ncatt_put(ncout, 0, "geospatial_lat_min", -90)
  ncatt_put(ncout, 0, "geospatial_lat_max",  90)
  ncatt_put(ncout, 0, "geospatial_lon_min", -180)
  ncatt_put(ncout, 0, "geospatial_lon_max",  180)
  
  ncatt_put(ncout, 0, "geospatial_vertical_min", 0)
  ncatt_put(ncout, 0, "geospatial_vertical_max", 0)
  
  ncatt_put(
    ncout,
    0,
    "time_coverage_start",
    time_coverage_start_val
  )
  
  ncatt_put(
    ncout,
    0,
    "time_coverage_end",
    time_coverage_end_val
  )
  
  ncatt_put(
    ncout,
    0,
    "time_coverage_duration",
    "P1M"
  )
  
  ncatt_put(
    ncout,
    0,
    "time_coverage_resolution",
    "P1M"
  )
  
  ncatt_put(
    ncout,
    0,
    "standard_name_vocabulary",
    "NetCDF Climate and Forecast (CF) Metadata Convention"
  )
  
  ncatt_put(
    ncout,
    0,
    "license",
    "ESA CCI Data Policy: free and open access"
  )
  
  ncatt_put(
    ncout,
    0,
    "platform",
    "Derived from existing burned area datasets based on Terra, Aqua, Sentinel-3 and Suomi-NPP"
  )
  
  ncatt_put(
    ncout,
    0,
    "sensor",
    "Derived from existing burned area datasets based on MODIS, VIIRS, OLCI and SLSTR"
  )
  
  ncatt_put(
    ncout,
    0,
    "spatial_resolution",
    "0.25 degrees"
  )
  
  ncatt_put(
    ncout,
    0,
    "key_variables",
    "burned_area"
  )
  
  ncatt_put(
    ncout,
    0,
    "geospatial_lon_units",
    "degrees_east"
  )
  
  ncatt_put(
    ncout,
    0,
    "geospatial_lat_units",
    "degrees_north"
  )
  
  ncatt_put(
    ncout,
    0,
    "geospatial_lon_resolution",
    0.25
  )
  
  ncatt_put(
    ncout,
    0,
    "geospatial_lat_resolution",
    0.25
  )
  
  # ==========================================================
  # 5.10) Escribir bounds
  # ==========================================================
  
  lon_b <- rbind(lon - half, lon + half)
  lat_b <- rbind(lat - half, lat + half)
  
  t_bnds <- as.numeric(
    c(tb[1], tb[2]) - as.Date("1970-01-01")
  )
  
  ncvar_put(
    nc = ncout,
    varid = "lon_bounds",
    vals = lon_b
  )
  
  ncvar_put(
    nc = ncout,
    varid = "lat_bounds",
    vals = lat_b
  )
  
  ncvar_put(
    nc = ncout,
    varid = "time_bounds",
    vals = array(t_bnds, dim = c(2, 1))
  )
  
  # ==========================================================
  # 5.11) Escribir datos
  # ==========================================================
  
  ba <- burned_area[, , k]
  ba[is.na(ba)] <- 0
  
  unc <- uncertainty[, , k]
  unc[is.na(unc)] <- 0
  
  ncvar_put(
    nc = ncout,
    varid = "burned_area",
    vals = ba,
    start = c(1, 1, 1),
    count = c(-1, -1, 1)
  )
  
  ncvar_put(
    nc = ncout,
    varid = "uncertainty",
    vals = unc,
    start = c(1, 1, 1),
    count = c(-1, -1, 1)
  )
  
  nc_close(ncout)
  
  message("✔ Escrito: ", out_path)
}

# ==========================================================
# 6) Ejecutar para todos los meses
# ==========================================================

for (k in seq_along(dates)) {
  write_month_nc(k)
}
# Paquetes
# install.packages("ncdf4")
# install.packages("fields")
library(ncdf4)
library(fields)
# Ruta del archivo (ej.: enero 2019)
f <- "/mnt/disco6tb/MRBA60-2/results/C1-MRBA60-2003-2024-NetCDF/2019/20190801-ESACCI-L4_FIRE-BA-MR_HAR-fv6.0.0.nc"
nc  <- nc_open(f)
lon <- ncvar_get(nc, "lon")
lat <- ncvar_get(nc, "lat")
ba  <- ncvar_get(nc, "burned_area")  # dims: lon x lat x time (1) → matriz [length(lon), length(lat)]
ba_units <- tryCatch(ncatt_get(nc, "burned_area", "units")$value, error = function(e) NA)
fillval  <- tryCatch(ncatt_get(nc, "burned_area", "_FillValue")$value, error = function(e) NA)
nc_close(nc)
rm(f)
ba=ba/1e6

ba[ba==0]=NA
# Limpia _FillValue a NA si aplica
if (is.numeric(fillval) && is.finite(fillval)) ba[ba == fillval] <- NA_real_

# Si 'ba' tiene dimensión extra de tiempo, reduce:
if (length(dim(ba)) == 3 && dim(ba)[3] == 1) ba <- ba[,,8]

# Nota: cuando pasas x=lon, y=lat, z debe ser matriz [length(lon) x length(lat)], ¡no transpongas!
par(mar = c(4, 4, 2, 5) + 0.1, xaxs = "i", yaxs = "i")
image.plot(
  x = lon, y = lat, z = (ba),
  xlab = "Longitud", ylab = "Latitud",
  main = "Burned Area — 2003-01 (FireCCI60)",
  legend.lab = ifelse(is.na(ba_units), "burned_area", paste0("burned_area (", ba_units, ")")),
  col = tim.colors(120),
  useRaster = TRUE
)
# ibrary(maps)

maps::map(
  "world",
  add = TRUE,
  col = "black",
  lwd = 0.6
)
# === Añadir bordes continentales ===
# map("world", add = TRUE, col = "black", lwd = 0.6)        # costas y fronteras simples
 # map("worldHires", package="mapdata", add = TRUE, col="black", lwd=0.5)  # más detalle

# dev.off()




