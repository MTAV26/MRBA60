rm(list = ls())
gc()

source("/mnt/disco6tb/MRBA60/scripts/0_functions/functions_fire_extract.R")

generar_rdata_fire(
  ruta_firecci51   = "/mnt/disco6tb/MRBA60/data/A1_RAW/FireCCI51_025degree-download/PSD_Grid",
  ruta_fireccis311 = "/mnt/disco6tb/MRBA60/data/A1_RAW/FireCCIS311_025degree-download/PSD_Grid",
  dir_out          = "/mnt/disco6tb/MRBA60/data/A3_ADJ",
  años_f51         = 2003:2024,
  años_s3          = 2019:2024
)

dir_adj <- "/mnt/disco6tb/MRBA60/data/A3_ADJ"

load(file.path(dir_adj, "longitude.RData"))
load(file.path(dir_adj, "latitude.RData"))
load(file.path(dir_adj, "FireCCI51_2003_2024_0.25degree.RData"))
load(file.path(dir_adj, "FireCCIS311_2019_2024_0.25degree.RData"))

print(dim(f51))
print(dim(s3))
print(length(lon))
print(length(lat))

factor_unidad <- 1e6
unidad_y <- "Burned Area (/1e6 m2)"

if (length(dim(f51)) != 3) stop("f51 no tiene 3 dimensiones.")
if (length(dim(s3))  != 3) stop("s3 no tiene 3 dimensiones.")

# FireCCI51 completo: 2003-2024
dates_f51 <- seq(as.Date("2003-01-01"), by = "month", length.out = dim(f51)[3])
# FireCCIS311 completo: 2019-2024
dates_s3  <- seq(as.Date("2019-01-01"), by = "month", length.out = dim(s3)[3])

ts_f51 <- apply(f51, 3, function(x) sum(x, na.rm = TRUE)) / factor_unidad
ts_s3  <- apply(s3,  3, function(x) sum(x, na.rm = TRUE)) / factor_unidad

df_f51 <- data.frame(date = dates_f51, f51 = ts_f51)
df_s3  <- data.frame(date = dates_s3,  s3  = ts_s3)

ylim_full <- range(c(df_f51$f51, df_s3$s3), na.rm = TRUE)

plot(
  df_f51$date, df_f51$f51,
  type = "l", lwd = 2, col = "black",
  ylim = ylim_full,
  xlab = "", ylab = unidad_y,
  main = "Monthly Total Burned Area"
)
lines(df_s3$date, df_s3$s3, col = "red", lwd = 2)

legend(
  "topright",
  legend = c("FireCCI51 (2003-2024)", "FireCCIS311 (2019-2024)"),
  col = c("black", "red"),
  lwd = 2,
  bty = "n"
)

# Periodo común actualizado: 2019-2024
idx_f51_common <- df_f51$date >= as.Date("2019-01-01") & df_f51$date <= as.Date("2024-12-01")
idx_s3_common  <- df_s3$date  >= as.Date("2019-01-01") & df_s3$date  <= as.Date("2024-12-01")

df_f51_common <- df_f51[idx_f51_common, ]
df_s3_common  <- df_s3[idx_s3_common, ]

ylim_common <- range(c(df_f51_common$f51, df_s3_common$s3), na.rm = TRUE)

plot(
  df_f51_common$date, df_f51_common$f51,
  type = "l", lwd = 2, col = "black",
  ylim = ylim_common,
  xlab = "", ylab = unidad_y,
  main = "Monthly Total Burned Area (2019-2024)"
)
lines(df_s3_common$date, df_s3_common$s3, col = "red", lwd = 2)

legend(
  "topright",
  legend = c("FireCCI51", "FireCCIS311"),
  col = c("black", "red"),
  lwd = 2,
  bty = "n"
)
