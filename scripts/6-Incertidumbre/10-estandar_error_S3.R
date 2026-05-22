# ============================================================
# Read Sentinel-3 SYN FireCCIS311 Standard Error (2019–2022)
# Build monthly array [lon, lat, time] with 48 months
# ============================================================

suppressPackageStartupMessages({
  library(ncdf4)
  library(lubridate)
})

# --- Directorios por año ---
dirs_years <- c(
  "2019" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2019",
  "2020" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2020",
  "2021" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2021",
  "2022" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2022",
  "2023" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2023",
  "2024" = "/home/miguel/discoX/0_Final_Products/FireCCIS311/PSD_Grid/2024"
)

# # --- Secuencia de fechas y archivos esperados ---
# dates_s3 <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")
# years_vec <- year(dates_s3)
# months_vec <- month(dates_s3)
# 
# files_expected <- character(length(dates_s3))
# for (i in seq_along(dates_s3)) {
#   y <- sprintf("%04d", years_vec[i])
#   m <- sprintf("%02d", months_vec[i])
#   d <- "01"
#   fname <- sprintf("%s%s%s-ESACCI-L4_FIRE-BA-SYN-fv1.1.nc", y, m, d)
#   files_expected[i] <- file.path(dirs_years[y], fname)
# }
# --- Secuencia de fechas y archivos esperados ---
dates_s3 <- seq(as.Date("2019-01-01"), as.Date("2024-12-01"), by = "month")

years_vec  <- year(dates_s3)
months_vec <- month(dates_s3)

files_expected <- character(length(dates_s3))

for (i in seq_along(dates_s3)) {
  
  y <- sprintf("%04d", years_vec[i])
  m <- sprintf("%02d", months_vec[i])
  d <- "01"
  
  # Julio y agosto de 2022 usan versión fv1.2.nc
  version <- if (years_vec[i] == 2022 && months_vec[i] %in% c(7, 8)) {
    "fv1.2internal"
  } else {
    "fv1.1"
  }
  
  fname <- sprintf(
    "%s%s%s-ESACCI-L4_FIRE-BA-SYN-%s.nc",
    y, m, d, version
  )
  
  files_expected[i] <- file.path(dirs_years[y], fname)
}
# --- Comprobar existencia de archivos ---
missing_files <- files_expected[!file.exists(files_expected)]
if (length(missing_files) > 0) {
  cat("⚠️ WARNING: faltan", length(missing_files), "archivos. Se omiten.\n")
  print(missing_files)
}

# --- Abrir el primer NetCDF disponible para dimensionar ---
first_idx <- which(file.exists(files_expected))[1]
if (is.na(first_idx)) stop("No se encontró ningún NetCDF de entrada.")

nc0 <- nc_open(files_expected[first_idx])

# Detectar variable de "standard error"
var_names <- names(nc0$var)
cand_idx <- grepl("standard.*error|std.*err|^se$|_se$|^rmse$|_rmse$", var_names, ignore.case = TRUE)
if (!any(cand_idx)) {
  cat("Variables en el archivo:\n"); print(var_names)
  stop("No se encontró una variable que parezca 'standard error'.")
}
# Elegimos la primera candidata (si hay varias, puedes ajustar aquí)
varname_se <- var_names[which(cand_idx)[1]]
cat("✅ Variable seleccionada para standard error:", varname_se, "\n")

# Leer lon/lat (intenta nombres comunes)
coord_lon_names <- c("lon","longitude","x")
coord_lat_names <- c("lat","latitude","y")
lon_name <- names(nc0$dim)[match(coord_lon_names, tolower(names(nc0$dim)) , nomatch = 0)][1]
lat_name <- names(nc0$dim)[match(coord_lat_names, tolower(names(nc0$dim)) , nomatch = 0)][1]
if (is.na(lon_name)) lon_name <- coord_lon_names[coord_lon_names %in% names(nc0$var)][1]
if (is.na(lat_name)) lat_name <- coord_lat_names[coord_lat_names %in% names(nc0$var)][1]

# Si estaban como variables:
if (!is.na(lon_name) && lon_name %in% names(nc0$var)) {
  lon_vec <- ncvar_get(nc0, lon_name)
} else {
  # si están como dimensiones:
  dim_lon <- nc0$dim[[ which(tolower(names(nc0$dim)) %in% coord_lon_names)[1] ]]
  lon_vec <- dim_lon$vals
}
if (!is.na(lat_name) && lat_name %in% names(nc0$var)) {
  lat_vec <- ncvar_get(nc0, lat_name)
} else {
  dim_lat <- nc0$dim[[ which(tolower(names(nc0$dim)) %in% coord_lat_names)[1] ]]
  lat_vec <- dim_lat$vals
}

# Dimensiones de la variable
vdims <- sapply(nc0$var[[varname_se]]$dim, function(d) d$name)
nlon <- length(lon_vec)
nlat <- length(lat_vec)

# Prealocar array [lon, lat, time]
SE_S3 <- array(NA_real_, dim = c(nlon, nlat, length(dates_s3)))

# Missing values
fill_att <- ncatt_get(nc0, varname_se, "_FillValue")$value
miss_att <- ncatt_get(nc0, varname_se, "missing_value")$value
nc_close(nc0)

# --- Lectura secuencial por mes ---
tpos <- 0L
for (i in seq_along(files_expected)) {
  f <- files_expected[i]
  tpos <- tpos + 1L
  if (!file.exists(f)) {
    cat("  · Saltando (no existe):", f, "\n")
    next
  }
  nc <- nc_open(f)
  v <- ncvar_get(nc, varname_se)
  # Convertir a NA los fill/missing conocidos
  if (!is.null(fill_att) && is.finite(fill_att)) v[v == fill_att] <- NA_real_
  if (!is.null(miss_att) && is.finite(miss_att)) v[v == miss_att] <- NA_real_
  
  # Alinear orientación a [lon,lat]
  dims_now <- sapply(nc$var[[varname_se]]$dim, function(d) d$name)
  # Si viene [lat,lon], transponer
  if (length(dim(v)) == 2) {
    if (tolower(dims_now[1]) %in% c("lat","latitude") && tolower(dims_now[2]) %in% c("lon","longitude")) {
      v <- t(v)  # pasa a [lon,lat]
    }
  } else if (length(dim(v)) > 2) {
    stop("La variable tiene más de 2 dimensiones; ajusta este bloque según el NetCDF.")
  }
  
  # Chequeo de tamaño
  if (!all(dim(v)[1] == nlon, dim(v)[2] == nlat)) {
    stop("Dimensiones inesperadas en ", basename(f), ": ", paste(dim(v), collapse=" x "),
         " (esperado: ", nlon, " x ", nlat, ")")
  }
  
  SE_S3[,,tpos] <- v
  nc_close(nc)
  if (i %% 6 == 0) { cat("  · Leídos", i, "de", length(files_expected), "archivos\n"); gc() }
}

# --- Dimensión temporal / nombres
dimnames(SE_S3) <- list(NULL, NULL, format(dates_s3, "%Y-%m"))

cat("✅ Construido SE_S3 con dimensiones:", paste(dim(SE_S3), collapse=" x "), "\n")

# Comprobación: ¿la latitud está decreciente?
if (length(lat_vec) < 2) stop("lat_vec tiene longitud < 2.")
is_decreasing <- lat_vec[2] < lat_vec[1]

if (is_decreasing) {
  cat("Latitud decreciente detectada. Reorientando latitud (sur->norte)...\n")
  lat_order <- rev(seq_along(lat_vec))   # índices invertidos
  lat_vec   <- lat_vec[lat_order]        # actualizar vector de latitudes
  SE_S3     <- SE_S3[, lat_order, , drop = FALSE]  # reordenar dimensión lat en el array
  cat("Hecho. Nueva latitud: de", lat_vec[1], "a", lat_vec[length(lat_vec)], "\n")
} else {
  cat("La latitud ya está en orden sur->norte. No se hacen cambios.\n")
}

library(maps)
library(fields)
# (Opcional) vista rápida para verificar orientación
image.plot(lon_vec, lat_vec, SE_S3[,,2], main = "SE_S3 (mes 1) - lat reorientada")



# output_dir_RData="C:/Users/migue/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-2/RData/"
output_dir_RData="/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/"

# --- Guardar a disco (ajusta ruta de salida si quieres otro sitio) ---
out_rdata <- file.path(output_dir_RData, "FireCCIS311_S3_SE_monthly_2019_2024.RData")
save(SE_S3, lon_vec, lat_vec, dates_s3, file = out_rdata)
cat("💾 Guardado:", out_rdata, "\n")

