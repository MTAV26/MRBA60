

rm(list = ls())
graphics.off()
gc()

output_dir_RData="/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/RData/"

#MASK TRUE/fALSE
# load(paste0(output_dir_RData, "MASK_NOHARMONISED.RData"))
# DATOS BA
load(paste0(output_dir_RData, "BA_Incertidumbre_MRBA60.RData"))
load(paste0(output_dir_RData, "FireCCIS311_S3_SE_monthly_2019_2024.RData"))


load("/mnt/disco6tb/MRBA60/data/A3_ADJ/longitude.RData")
load("/mnt/disco6tb/MRBA60/data/A3_ADJ/latitude.RData")
dim(incert_Fire60)
dim(SE_S3)

summary(as.vector(incert_Fire60))
BA_60_192<-incert_Fire60[,,1:192]
summary(as.vector(BA_60_192))

# BA_60_192[is.na(BA_60_192)]=0
# rm(incert_Fire60);gc()
BA_60_192=BA_60_192*1e6

# stopifnot(identical(dim(BA_60_192)[1:2], dim(BA_S3)[1:2]))

nxy <- dim(BA_60_192)[1:2]
out <- array(NA_real_, dim = c(nxy, 192 + 72))  # 1440 x 720 x 240
dim(out)
out[,, 1:192]   <- BA_60_192
out[,, 193:264] <- SE_S3

summary(as.vector(SE_S3))
summary(as.vector(BA_60_192))

image.plot(BA_60_192[,,1])
image.plot(SE_S3[,,1])


Unc_MRBA60H=out
# # --- Guardar a disco (ajusta ruta de salida si quieres otro sitio) ---
# out_rdata <- file.path(output_dir_RData, "FireCCII60_Unc_m2_monthly_2003_2024.RData")
# save(Unc_FireCCI60, lon, lat, file = out_rdata)
# cat("💾 Guardado:", out_rdata, "\n")
# 
# 
# 
# load(paste0(output_dir_RData, "FireCCII60_Unc_m2_monthly_2003_2024.RData"))
# summary(as.vector(Unc_FireCCI60))
# 






Unc_MRBA60H=Unc_MRBA60H/1e6
# ---------------------------
# MÁSCARAS FÍSICAS (lat-dependientes)
# ---------------------------
ncol_grid <- length(lon); nrow_grid <- length(lat)
cell_area_constant <- (110.57 * 0.25) * (111.32 * 0.25)  # ≈ 769.29 km² en el ecuador
area_by_row <- cell_area_constant * cos(lat * pi/180)
area_matrix <- t(matrix(area_by_row, ncol = ncol_grid, nrow = nrow_grid, byrow = FALSE))
image.plot(lon, lat, area_matrix)
# amin_px_eq <- 0.09
# amin_by_row <- amin_px_eq * cos(lat * pi/180)
# amin_matrix <- t(matrix(amin_by_row, ncol = ncol_grid, nrow = nrow_grid, byrow = FALSE))
dev.off()

ncols <- length(lon)   # debería ser 1440
nrows <- length(lat)   # debería ser 720
amin_px_eq  <- 0.09
amin_by_row <- amin_px_eq * cos(lat * pi/180)        # length = nrows
amin_matrix <- t(matrix(amin_by_row, nrow = nrows, ncol = ncols, byrow = FALSE))
image.plot(lon, lat, amin_matrix)
dev.off()

# Unc_FireCCI60: [lon, lat, time]  (1440 x 720 x 264)
# area_matrix  : [lon, lat]        (1440 x 720)

# Replicar area_matrix a 3D (mismo tamaño que Unc_FireCCI60)
area_cap <- array(area_matrix, dim = dim(Unc_MRBA60H))
# Máscara de celdas donde la incertidumbre supera el área máxima (y ambos no son NA)
mask_cap <- !is.na(Unc_MRBA60H) & !is.na(area_cap) & (Unc_MRBA60H > area_cap)
# Aplicar el tope
Unc_MRBA60H[mask_cap] <- area_cap[mask_cap]
# (Opcional) resumen de cambios
cat("Celdas/meses capadas:", sum(mask_cap), "\n")
summary(as.vector(Unc_MRBA60H))
dim(amin_matrix)

# Unc_MRBA60H: [lon, lat, time]  (1440 x 720 x 264)
# amin_matrix  : [lon, lat]        (1440 x 720)
# Replicar amin_matrix a 3D
amin_cap <- array(amin_matrix, dim = dim(Unc_MRBA60H))
# Máscara: incertidumbre por debajo del mínimo permitido (y sin NA en el umbral)
mask_low <- !is.na(Unc_MRBA60H) & !is.na(amin_cap) & (Unc_MRBA60H < amin_cap)
# Forzar a 0 donde esté por debajo del mínimo
Unc_MRBA60H[mask_low] <- 0
# (Opcional) resumen
cat("Celdas/meses forzadas a 0 por estar < mínimo:", sum(mask_low), "\n")
Unc_MRBA60=Unc_MRBA60H*1e6

# --- Guardar a disco (ajusta ruta de salida si quieres otro sitio) ---
out_rdata <- file.path(output_dir_RData, "MRBA60_Unc_m2_monthly_2003_2024.RData")
save(Unc_MRBA60, lon, lat, file = out_rdata)
cat("💾 Guardado:", out_rdata, "\n")


# summary(as.vector(Unc_FireCCI60[,,1]))

# area_matrix=area_matrix*1e6

# summary(as.vector(area_matrix))
