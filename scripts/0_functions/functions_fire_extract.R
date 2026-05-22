rm(list = ls())
gc()

library(terra)

# =========================
# FUNCIÓN: extraer burned_area de NetCDF
# =========================
extraer_burned_area <- function(stack_obj, nombre_obj = "objeto") {
  ba_indices <- grep("^burned_area$", names(stack_obj))
  
  if (length(ba_indices) == 0) {
    ba_indices <- grep("burned_area", names(stack_obj), ignore.case = TRUE)
  }
  
  if (length(ba_indices) == 0) {
    stop("No se encontraron capas 'burned_area' en ", nombre_obj)
  }
  
  ba_stack <- stack_obj[[ba_indices]]
  return(ba_stack)
}

# =========================
# FUNCIÓN: leer un archivo mensual
# Admite .nc y .tif
# =========================
leer_archivo_mensual_fire <- function(archivo) {
  ext <- tolower(tools::file_ext(archivo))
  r <- rast(archivo)
  
  if (ext == "nc") {
    r <- extraer_burned_area(r, basename(archivo))
  } else if (ext == "tif") {
    if (!grepl("_m2\\.tif$", basename(archivo), ignore.case = TRUE)) {
      stop(
        "El archivo TIFF no parece estar en m2: ", basename(archivo),
        ". Se esperaba un nombre terminado en '_m2.tif'."
      )
    }
  } else {
    stop("Formato no soportado: ", archivo)
  }
  
  if (nlyr(r) != 1) {
    stop(
      "El archivo mensual ", basename(archivo),
      " debe aportar una única capa tras la lectura, pero tiene ", nlyr(r), "."
    )
  }
  
  names(r) <- tools::file_path_sans_ext(basename(archivo))
  return(r)
}

# =========================
# FUNCIÓN: listar archivos por año
# =========================
listar_archivos_fire <- function(ruta, años = NULL) {
  carpetas <- list.dirs(ruta, recursive = FALSE, full.names = TRUE)
  
  if (length(carpetas) == 0) {
    stop("No se encontraron subcarpetas en la ruta: ", ruta)
  }
  
  if (!is.null(años)) {
    carpetas <- carpetas[basename(carpetas) %in% as.character(años)]
  }
  
  if (length(carpetas) == 0) {
    stop("No se encontraron carpetas para los años solicitados en: ", ruta)
  }
  
  carpetas <- carpetas[order(as.integer(basename(carpetas)))]
  
  archivos <- unlist(lapply(carpetas, function(carpeta) {
    f <- list.files(
      carpeta,
      pattern = "\\.(nc|tif)$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    f <- sort(f)
    
    if (length(f) != 12) {
      stop(
        "La carpeta ", carpeta,
        " no tiene 12 archivos mensuales. Encontrados: ", length(f)
      )
    }
    
    return(f)
  }))
  
  if (length(archivos) == 0) {
    stop("No se encontraron archivos .nc o .tif en la ruta: ", ruta)
  }
  
  return(archivos)
}

# =========================
# FUNCIÓN: leer serie completa y validar geometría
# =========================
procesar_firecci <- function(ruta, años = NULL, nombre_producto = "producto") {
  archivos <- listar_archivos_fire(ruta, años = años)
  
  cat("Número de archivos encontrados para ", nombre_producto, ": ", length(archivos), "\n", sep = "")
  
  rasters <- vector("list", length(archivos))
  rasters[[1]] <- leer_archivo_mensual_fire(archivos[1])
  ref <- rasters[[1]]
  
  for (i in seq_along(archivos)) {
    if (i > 1) {
      rasters[[i]] <- leer_archivo_mensual_fire(archivos[i])
      
      if (!compareGeom(ref, rasters[[i]], stopOnError = FALSE)) {
        stop(
          "Geometría distinta en ", nombre_producto,
          ": ", basename(archivos[1]), " vs ", basename(archivos[i])
        )
      }
    }
  }
  
  stack <- do.call(c, rasters)
  return(stack)
}

# =========================
# FUNCIÓN: reorganizar array
# =========================
array_transpuesta <- function(mat) {
  arr <- aperm(mat, c(2, 1, 3))
  arr <- arr[, ncol(arr):1, ]
  return(arr)
}

# =========================
# FUNCIÓN: generar y guardar RData
# =========================
generar_rdata_fire <- function(
    ruta_firecci51 = "/mnt/disco6tb/MHBA60/data/A1_RAW/FireCCI51_025degree-download/PSD_Grid",
    ruta_fireccis311 = "/mnt/disco6tb/MHBA60/data/A1_RAW/FireCCIS311_025degree-download/PSD_Grid",
    dir_out = "/mnt/disco6tb/MHBA60/data/A3_ADJ",
    años_f51 = 2003:2024,
    años_s3 = 2019:2024
) {
  dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)
  
  cat("Leyendo FireCCI51...\n")
  stack_f51 <- procesar_firecci(
    ruta_firecci51,
    años = años_f51,
    nombre_producto = "FireCCI51"
  )
  
  cat("Leyendo FireCCIS311...\n")
  stack_s3 <- procesar_firecci(
    ruta_fireccis311,
    años = años_s3,
    nombre_producto = "FireCCIS311"
  )
  
  if (!compareGeom(stack_f51[[1]], stack_s3[[1]], stopOnError = FALSE)) {
    stop("FireCCI51 y FireCCIS311 no tienen la misma geometría base.")
  }
  
  ext_obj <- ext(stack_f51)
  res_x <- res(stack_f51)[1]
  res_y <- res(stack_f51)[2]
  
  lon <- seq(ext_obj[1] + res_x / 2, ext_obj[2] - res_x / 2, by = res_x)
  lat <- seq(ext_obj[3] + res_y / 2, ext_obj[4] - res_y / 2, by = res_y)
  
  cat("Extrayendo burned area / capas mensuales de FireCCI51...\n")
  f51 <- array_transpuesta(as.array(stack_f51))
  
  cat("Extrayendo burned area / capas mensuales de FireCCIS311...\n")
  s3 <- array_transpuesta(as.array(stack_s3))
  
  cat("Comprobando dimensiones...\n")
  if (!identical(dim(f51)[1:2], dim(s3)[1:2])) {
    stop("FireCCI51 y FireCCIS311 no coinciden en dimensión espacial.")
  }
  
  cat("Guardando RData...\n")
  save(lon, file = file.path(dir_out, "longitude.RData"))
  save(lat, file = file.path(dir_out, "latitude.RData"))
  save(f51, file = file.path(dir_out, "FireCCI51_2003_2024_0.25degree.RData"))
  save(s3,  file = file.path(dir_out, "FireCCIS311_2019_2024_0.25degree.RData"))
  
  cat("Dimensiones f51:", paste(dim(f51), collapse = " x "), "\n")
  cat("Dimensiones s3 :", paste(dim(s3), collapse = " x "), "\n")
  cat("Longitud lon   :", length(lon), "\n")
  cat("Longitud lat   :", length(lat), "\n")
  cat("Capas f51      :", dim(f51)[3], "\n")
  cat("Capas s3       :", dim(s3)[3], "\n")
  
  rm(stack_f51, stack_s3, ext_obj, res_x, res_y, lon, lat, f51, s3)
  gc()
  
  cat("Proceso finalizado.\n")
}
