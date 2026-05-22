# ===============================
# CONFIGURACIÓN PREVIA
# ===============================

rm(list = ls())
graphics.off()
gc()

library(ncdf4)

# ============================================================================
# DIRECTORIOS
# ============================================================================

dir_adj <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"
dir_tem <- "/mnt/disco6tb/MRBA60/data/A2_TEMP/"
dir_mask_fire <- "/mnt/disco6tb/MRBA60/results/A1-Built-mask/MRBA60-2003-2024-V1/MASK_FIRE_PREPROC/RData/"

# ============================================================================
# ARCHIVOS DE ENTRADA
# ============================================================================

file_landsea <- paste0(dir_tem, "land_sea_mask_025degree_binary_1440x720.nc")
file_firemask <- paste0(dir_mask_fire, "FireMask_AF3030F.RData")

# ============================================================================
# ARCHIVOS DE SALIDA
# ============================================================================

file_landsea_adj_nc <- paste0(
  dir_adj,
  "land_sea_mask_025degree_binary_1440x720_fire_adjusted.nc"
)

file_landsea_adj_rdata <- paste0(
  dir_adj,
  "land_sea_mask_025degree_binary_1440x720_fire_adjusted.RData"
)

# ============================================================================
# CARGAR LAND/SEA MASK
# ============================================================================

nc_land <- nc_open(file_landsea)

print(names(nc_land$var))

lon_land <- ncvar_get(nc_land, "lon")
lat_land <- ncvar_get(nc_land, "lat")

land_sea_mask <- ncvar_get(nc_land, "sftlf")

nc_close(nc_land)

cat("Dimensiones land_sea_mask original:\n")
print(dim(land_sea_mask))

cat("Rango lon_land:\n")
print(range(lon_land, na.rm = TRUE))

cat("Rango lat_land:\n")
print(range(lat_land, na.rm = TRUE))

cat("Valores land_sea_mask original:\n")
print(table(as.vector(land_sea_mask), useNA = "ifany"))

# ============================================================================
# CARGAR FIRE MASK
# ============================================================================

load(file_firemask)

mask_fire <- FireMask_AF3030F

cat("\nDimensiones mask_fire:\n")
print(dim(mask_fire))

cat("Valores mask_fire:\n")
print(table(as.vector(mask_fire), useNA = "ifany"))


# ============================================================================
# IDENTIFICAR CELDAS CON FUEGO ALGUNA VEZ
# ============================================================================

# TRUE si una celda tuvo fuego al menos una vez en 2003-2024
fire_ever <- apply(mask_fire == 1, c(1, 2), any, na.rm = TRUE)

cat("\nDimensiones fire_ever:\n")
print(dim(fire_ever))

n_fire_ever <- sum(fire_ever, na.rm = TRUE)

cat("Número de celdas espaciales con fuego alguna vez:\n")
print(n_fire_ever)

# ============================================================================
# IDENTIFICAR CELDAS DE OCÉANO QUE TUVIERON FUEGO
# ============================================================================

# Convención:
# land_sea_mask == 0 -> océano
# land_sea_mask == 1 -> tierra

ocean_with_fire <- fire_ever & land_sea_mask == 0

n_ocean_with_fire <- sum(ocean_with_fire, na.rm = TRUE)

# Recuentos base
n_ocean_total <- sum(land_sea_mask == 0, na.rm = TRUE)
n_land_total  <- sum(land_sea_mask == 1, na.rm = TRUE)
n_grid_total  <- sum(!is.na(land_sea_mask))

# Porcentajes
perc_ocean_adapted <- 100 * n_ocean_with_fire / n_ocean_total
perc_grid_adapted  <- 100 * n_ocean_with_fire / n_grid_total
perc_fire_ever_ocean <- 100 * n_ocean_with_fire / n_fire_ever

cat("\n================ RESUMEN DE CONFLICTOS ================\n")

cat("Celdas clasificadas como océano pero con fuego alguna vez:\n")
print(n_ocean_with_fire)

cat("\nTotal celdas océano en land_sea_mask original:\n")
print(n_ocean_total)

cat("\nTotal celdas tierra en land_sea_mask original:\n")
print(n_land_total)

cat("\nTotal celdas válidas del grid:\n")
print(n_grid_total)

cat("\nTotal celdas con fuego alguna vez:\n")
print(n_fire_ever)

cat("\nPorcentaje de celdas de océano reclasificadas a tierra (%):\n")
print(round(perc_ocean_adapted, 6))

cat("\nPorcentaje respecto al total de celdas del grid (%):\n")
print(round(perc_grid_adapted, 6))

cat("\nPorcentaje de celdas con fuego alguna vez que estaban clasificadas como océano (%):\n")
print(round(perc_fire_ever_ocean, 6))

# ============================================================================
# MAPA RÁPIDO DE CONFLICTOS
# ============================================================================

image(
  ocean_with_fire,
  main = "Celdas océano en land_sea con fuego alguna vez",
  xlab = "lon index",
  ylab = "lat index"
)

# ============================================================================
# AJUSTAR LAND/SEA MASK
# ============================================================================

land_sea_mask_adj <- land_sea_mask

# Donde land_sea decía océano pero FireMask tuvo fuego alguna vez,
# se cambia a tierra
land_sea_mask_adj[ocean_with_fire] <- 1

# Asegurar valores binarios enteros
land_sea_mask_adj <- round(land_sea_mask_adj)
land_sea_mask_adj[land_sea_mask_adj != 0 & land_sea_mask_adj != 1] <- NA

cat("\n================ AJUSTE LAND/SEA ================\n")

cat("Valores land_sea_mask original:\n")
print(table(as.vector(land_sea_mask), useNA = "ifany"))

cat("\nValores land_sea_mask ajustada:\n")
print(table(as.vector(land_sea_mask_adj), useNA = "ifany"))

cat("\nNúmero de celdas cambiadas de océano a tierra:\n")
print(sum(land_sea_mask == 0 & land_sea_mask_adj == 1, na.rm = TRUE))

# ============================================================================
# COMPROBACIÓN FINAL
# ============================================================================

# Después del ajuste, no debería quedar ninguna celda con fuego alguna vez
# clasificada como océano

remaining_ocean_with_fire <- fire_ever & land_sea_mask_adj == 0

n_remaining_ocean_with_fire <- sum(remaining_ocean_with_fire, na.rm = TRUE)

cat("\n================ COMPROBACIÓN FINAL ================\n")

cat("Celdas con fuego alguna vez que siguen como océano tras el ajuste:\n")
print(n_remaining_ocean_with_fire)

if (n_remaining_ocean_with_fire == 0) {
  cat("OK: todas las celdas con fuego alguna vez están clasificadas como tierra.\n")
} else {
  warning("Todavía quedan celdas con fuego clasificadas como océano.")
}

# ============================================================================
# GUARDAR COMO RDATA
# ============================================================================

save(
  land_sea_mask_adj,
  lon_land,
  lat_land,
  fire_ever,
  ocean_with_fire,
  remaining_ocean_with_fire,
  n_ocean_with_fire,
  n_ocean_total,
  n_land_total,
  n_grid_total,
  n_fire_ever,
  n_remaining_ocean_with_fire,
  perc_ocean_adapted,
  perc_grid_adapted,
  perc_fire_ever_ocean,
  file = file_landsea_adj_rdata
)

cat("\nArchivo RData guardado en:\n")
cat(file_landsea_adj_rdata, "\n")

# ============================================================================
# GUARDAR COMO NETCDF
# ============================================================================

# Definir dimensiones
dim_lon <- ncdim_def(
  name = "lon",
  units = "degrees_east",
  vals = lon_land
)

dim_lat <- ncdim_def(
  name = "lat",
  units = "degrees_north",
  vals = lat_land
)

# Definir variable
var_land <- ncvar_def(
  name = "sftlf",
  units = "1",
  dim = list(dim_lon, dim_lat),
  missval = -9999,
  longname = "Adjusted binary land sea mask, ocean cells with fire detection set to land",
  prec = "short"
)

# Crear NetCDF
nc_out <- nc_create(file_landsea_adj_nc, vars = list(var_land))

# Escribir variable
ncvar_put(nc_out, var_land, land_sea_mask_adj)

# Atributos informativos de variable
ncatt_put(nc_out, "sftlf", "flag_values", c(0, 1), prec = "short")
ncatt_put(nc_out, "sftlf", "flag_meanings", "ocean land")
ncatt_put(
  nc_out,
  "sftlf",
  "description",
  "Binary land-sea mask adjusted using FireMask_AF3030F. Cells classified as ocean but with fire detected at least once were set to land."
)

# Atributos globales
ncatt_put(nc_out, 0, "source_land_sea_mask", basename(file_landsea))
ncatt_put(nc_out, 0, "source_fire_mask", basename(file_firemask))
ncatt_put(
  nc_out,
  0,
  "adjustment",
  "land_sea_mask set to 1 where FireMask_AF3030F == 1 at least once over time"
)

ncatt_put(nc_out, 0, "n_ocean_with_fire", n_ocean_with_fire)
ncatt_put(nc_out, 0, "n_ocean_total", n_ocean_total)
ncatt_put(nc_out, 0, "n_land_total", n_land_total)
ncatt_put(nc_out, 0, "n_grid_total", n_grid_total)
ncatt_put(nc_out, 0, "n_fire_ever", n_fire_ever)

ncatt_put(nc_out, 0, "perc_ocean_adapted", perc_ocean_adapted)
ncatt_put(nc_out, 0, "perc_grid_adapted", perc_grid_adapted)
ncatt_put(nc_out, 0, "perc_fire_ever_ocean", perc_fire_ever_ocean)

nc_close(nc_out)

cat("\nArchivo NetCDF guardado en:\n")
cat(file_landsea_adj_nc, "\n")


image(
  land_sea_mask_adj,
  main = "Celdas océano en land_sea con fuego alguna vez",
  xlab = "lon index",
  ylab = "lat index"
)

