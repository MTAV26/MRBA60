# ==== Paquetes ====
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


# ---------------------------
# Modelo <- "HBA-RF-2003-2022-CL60-R30"
Modelo <- "B1-MRBA60-2003-2024"

dir_oss <- "/mnt/disco6tb/MRBA60/data/A3_ADJ/"

output_dir       <- file.path("/mnt/disco6tb/MRBA60-2/results", Modelo)
output_dir_csv <- paste0(output_dir, "/csv/")
output_dir_plot <- paste0(output_dir, "/plot/")
# output_dir_plot_rle <- paste0(output_dir, "/plot_rle/")
output_dir_RData <- paste0(output_dir, "/RData/")



# ==== Configura tu carpeta de CSV ====
# Usa la que corresponda a tu modelo
dir_csv <- output_dir_csv #"/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/FireCCI60-2003-2022-CL30-R30-3/csv/"

# ==== Listar archivos ====
files <- list.files(dir_csv, pattern = "^SelectedPredictors_.*\\.csv$", full.names = TRUE)
stopifnot("No SelectedPredictors_*.csv files found" = length(files) > 0)

# ==== Diccionario de nombres bonitos ====
rename_map <- c(
  "f5"               = "FireCCI51",
  "count_ActiveFire" = "Active Fire",
  "prec"             = "Precipitation",
  "temp"             = "Temperature",
  "FRPsum"           = "FRP sum",
  "FRPmedian"        = "FRP median",
  "FWI"              = "FWId95",
  "lat"              = "Latitude",
  "lon"              = "Longitude",
  # "LAIHIGH"          = "LAI High",
  # "LAILOW"           = "LAI Low",
  "wind"             = "Wind Speed",
  "NDVI"             = "NDVI",
  "cloud"            = "Cloud",
  "gpp"              = "GPP",
  "vpd"              = "VPD",
  "soil"             = "Soil Moisture"
)

# ==== Lector robusto por archivo ====
read_predictors_from_csv <- function(fpath) {
  df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  
  # Detectar columnas (case-insensitive)
  month_col <- names(df)[grepl("^month$", names(df), ignore.case = TRUE)]
  # Prioriza RFE_*select* y excluye "Formula"
  rfe_candidates <- names(df)[grepl("^rfe", names(df), ignore.case = TRUE) &
                                !grepl("formula", names(df), ignore.case = TRUE)]
  rfe_sel <- rfe_candidates[grepl("select", rfe_candidates, ignore.case = TRUE)]
  rfe_col <- if (length(rfe_sel)) rfe_sel[1] else if (length(rfe_candidates)) rfe_candidates[1] else character(0)
  
  if (length(month_col) == 0 || length(rfe_col) == 0) {
    warning("Skipping file without Month/RFE: ", basename(fpath))
    return(NULL)
  }
  month_col <- month_col[1]
  
  # Nombre del bioma desde el archivo
  bio <- basename(fpath)
  bio <- sub("^SelectedPredictors_(.*)\\.csv$", "\\1", bio)
  bio <- gsub("__", " ", bio)
  bio <- gsub("_",  " ", bio)
  bio <- trimws(bio)
  
  # Extraer y limpiar
  df2 <- df %>%
    dplyr::transmute(
      Month = suppressWarnings(as.integer(trimws(.data[[month_col]]))),
      RFE   = as.character(.data[[rfe_col]])
    ) %>%
    dplyr::mutate(
      RFE   = ifelse(is.na(RFE), "", RFE),
      Biome = bio
    ) %>%
    dplyr::filter(!is.na(Month), Month >= 1, Month <= 12) %>%
    # Separar por + , o ; (con o sin espacios)
    tidyr::separate_rows(RFE, sep = "\\s*[+,;]\\s*") %>%
    dplyr::mutate(RFE = trimws(RFE)) %>%
    dplyr::filter(RFE != "")
  
  if (nrow(df2) == 0) return(NULL)
  df2
}

# ==== Leer todos los archivos ====
lst <- lapply(files, read_predictors_from_csv)
lst <- lst[!vapply(lst, is.null, logical(1))]
df_preds <- dplyr::bind_rows(lst)

# ==== Renombrar predictores (los no mapeados quedan igual) ====
df_preds <- df_preds %>%
  dplyr::mutate(RFE = dplyr::recode(RFE, !!!rename_map, .default = RFE))

# ==== Conteo por predictor y mes (nº de biomas donde entra) ====
conteo_preds <- df_preds %>%
  dplyr::group_by(Month, RFE) %>%
  dplyr::summarise(N_biomes = dplyr::n_distinct(Biome), .groups = "drop") %>%
  dplyr::mutate(
    Month_lab = factor(month.abb[Month], levels = month.abb),
    text_col  = ifelse(N_biomes <= 25, "black", "white")  # color del número
  )

# ==== (Opcional) limitar a los Top-N predictores por frecuencia total ====
# topN <- 30
# top_predictors <- conteo_preds %>%
#   dplyr::group_by(RFE) %>%
#   dplyr::summarise(Total = sum(N_biomes), .groups = "drop") %>%
#   dplyr::slice_max(Total, n = topN) %>%
#   dplyr::pull(RFE)
# conteo_preds_plot <- conteo_preds %>% dplyr::filter(RFE %in% top_predictors)
# (Si no quieres filtrar, usa 'conteo_preds' directamente:)
# conteo_preds_plot <- conteo_preds
# ==== Ordenar predictores por frecuencia total (menor -> mayor) ====
order_predictors <- conteo_preds %>%
  dplyr::group_by(RFE) %>%
  dplyr::summarise(Total = sum(N_biomes), .groups = "drop") %>%
  dplyr::arrange(Total) %>%
  dplyr::pull(RFE)

conteo_preds_plot <- conteo_preds %>%
  dplyr::mutate(RFE = factor(RFE, levels = order_predictors))

# ==== Tabla wide predictor × mes ====
sub <- as.data.frame(conteo_preds)[, c("RFE", "Month_lab", "N_biomes"), drop = FALSE]
conteo_wide <- tidyr::pivot_wider(
  sub,
  names_from  = Month_lab,
  values_from = N_biomes,
  values_fill = 0
)
conteo_wide <- conteo_wide[order(conteo_wide$RFE), ]
print(conteo_wide)
p_pred <- ggplot(conteo_preds_plot, aes(x = Month_lab, y = RFE, fill = N_biomes)) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = N_biomes, color = text_col), size = 5) +   # ← más grande
  scale_color_identity() +
  scale_fill_gradient(low = "grey90", high = "black", name = "") +
  labs(
    x = "", y = "Predictors (ordered by total frequency)",
    title = "Predictor selection frequency by month across biomes",
    subtitle = "Counts = number of biomes where predictor is selected"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black", face = "bold"),
    axis.title.y = element_text(size = 16, colour = "black", face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(2.5, "cm"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.text  = element_text(size = 12, colour = "black")
  )

print(p_pred)

print(p_pred)



# ==== (Opcional) exportar ====
# write.csv(conteo_wide, file.path(dir_csv, "predictor_month_counts_wide.csv"), row.names = FALSE)
output_dir= "/mnt/disco6tb/Dropbox/UAH/FireCCI60/ATBD/Figure6/"
ggsave(file.path(output_dir, "matrix_conteo_conteo_predictores.png"),
       plot = p_pred, width = 10, height = 12, dpi = 300)

# ==== Paquetes ====
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# ==== Directorio de trabajo ====
# # dir_csv <- "/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/M2_RF-2003-2022-CL60-R30/csv/"
# 
# # ==== Diccionario nombres bonitos ====
# rename_map <- c(
#   "f5"               = "FireCCI51",
#   "count_ActiveFire" = "Active Fire",
#   "prec"             = "Precipitation",
#   "temp"             = "Temperature",
#   "FRPsum"           = "FRP sum",
#   "FRPmedian"        = "FRP median",
#   "FWI"              = "FWId95",
#   "lat"              = "Latitude",
#   "LAIHIGH"          = "LAI High",
#   "LAILOW"           = "LAI Low",
#   "wind"             = "Wind Speed",
#   "NDVI"             = "NDVI",
#   "cloud"            = "Cloud",
#   "gpp"              = "GPP",
#   "vpd"              = "VPD",
#   "soil"             = "SOIL"
# )
# ==== Función lectora robusta para SHAP ====
read_shap_from_csv <- function(fpath) {
  
  df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
  
  # Limpiar nombres vacíos
  names(df) <- trimws(names(df))
  empty_names <- which(is.na(names(df)) | names(df) == "")
  if (length(empty_names) > 0) {
    names(df)[empty_names] <- paste0("X", empty_names)
  }
  
  # Detectar columnas
  month_col <- names(df)[grepl("^month$", names(df), ignore.case = TRUE)]
  var_col   <- names(df)[grepl("^variable$|feature|predictor", names(df), ignore.case = TRUE)]
  rank_col  <- names(df)[grepl("^rank$", names(df), ignore.case = TRUE)]
  
  if (length(month_col) == 0 || length(var_col) == 0 || length(rank_col) == 0) {
    warning("Saltando archivo sin Month/Variable/Rank: ", basename(fpath))
    return(NULL)
  }
  
  month_col <- month_col[1]
  var_col   <- var_col[1]
  rank_col  <- rank_col[1]
  
  # Nombre de bioma desde columna Biome si existe; si no, desde el archivo
  if ("Biome" %in% names(df)) {
    bio_vec <- as.character(df$Biome)
  } else {
    bio <- basename(fpath)
    bio <- sub("^shap_summary_(.*)\\.csv$", "\\1", bio, ignore.case = TRUE)
    bio <- sub("^SHAP_RF_(.*)_COMMON\\.csv$", "\\1", bio, ignore.case = TRUE)
    bio <- gsub("__", " ", bio)
    bio <- gsub("_",  " ", bio)
    bio <- trimws(bio)
    bio_vec <- bio
  }
  
  df2 <- df %>%
    dplyr::transmute(
      Month = suppressWarnings(as.integer(.data[[month_col]])),
      RFE   = as.character(.data[[var_col]]),
      Rank  = suppressWarnings(as.integer(.data[[rank_col]])),
      Biome = bio_vec
    ) %>%
    dplyr::filter(!is.na(Month), Month >= 1, Month <= 12)
  
  if (nrow(df2) == 0) return(NULL)
  
  df2
}

# ==== Leer todos los archivos SHAP ====
files_shap <- list.files(
  dir_csv,
  pattern = "^(shap_summary_|SHAP_RF_).*\\.csv$",
  full.names = TRUE,
  ignore.case = TRUE
)

if (length(files_shap) == 0) {
  stop("No se encontraron archivos SHAP en dir_csv: ", dir_csv)
}

lst_shap <- lapply(files_shap, read_shap_from_csv)
lst_shap <- lst_shap[!vapply(lst_shap, is.null, logical(1))]

if (length(lst_shap) == 0) {
  stop("Se encontraron archivos, pero ninguno contenía columnas Month, Variable y Rank.")
}

df_shap <- dplyr::bind_rows(lst_shap)

if (!"RFE" %in% names(df_shap)) {
  stop("La columna RFE no existe después de leer los archivos. Revisa nombres de columnas.")
}

# ==== Renombrar predictores ====
df_shap <- df_shap %>%
  dplyr::mutate(
    RFE = dplyr::recode(RFE, !!!rename_map, .default = RFE)
  )

# ==== Conteo: total y Rank = 1 ====
conteo_shap <- df_shap %>%
  dplyr::group_by(Month, RFE) %>%
  dplyr::summarise(
    N_total = dplyr::n_distinct(Biome),
    N_rank1 = sum(Rank == 1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Month_lab = factor(month.abb[Month], levels = month.abb),
    Label = paste0(N_total, " (", N_rank1, ")"),
    text_col = ifelse(N_total <= 25, "black", "white")
  )

# ==== Ordenar predictores por veces que fueron Rank 1 ====
order_predictors_rank <- conteo_shap %>%
  dplyr::group_by(RFE) %>%
  dplyr::summarise(Total_rank1 = sum(N_rank1), .groups = "drop") %>%
  dplyr::arrange(Total_rank1) %>%
  dplyr::pull(RFE)

conteo_shap <- conteo_shap %>%
  dplyr::mutate(
    RFE = factor(RFE, levels = order_predictors_rank)
  )

# ==== Heatmap combinado ====
p_shap <- ggplot(conteo_shap, aes(x = Month_lab, y = RFE, fill = N_total)) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = Label, color = text_col), size = 4) +
  scale_color_identity() +
  scale_fill_gradient(low = "grey90", high = "black", name = "") +
  labs(
    x = "Month",
    y = "Predictors (ordered by #1 importance)",
    title = "Predictor frequency and top-rank importance by month",
    subtitle = "Label = Total count (times ranked 1)"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black", face = "bold"),
    axis.title.y = element_text(size = 16, colour = "black", face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, colour = "black"),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(2.5, "cm"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.text = element_text(size = 12, colour = "black")
  )

print(p_shap)
# print(p_shap)
dev.off()


library(dplyr)
library(ggplot2)

# --- nº de biomas por mes (denominador) ---
biomes_per_month <- df_shap %>%
  dplyr::group_by(Month) %>%
  dplyr::summarise(Biomes_m = dplyr::n_distinct(Biome), .groups = "drop")

# --- conteo y porcentajes ---
conteo_shap <- df_shap %>%
  dplyr::group_by(Month, RFE) %>%
  dplyr::summarise(
    N_total = dplyr::n_distinct(Biome),     # cuántos biomas lo incluyen
    N_rank1 = sum(Rank == 1, na.rm = TRUE), # cuántos biomas fue #1
    .groups = "drop"
  ) %>%
  dplyr::left_join(biomes_per_month, by = "Month") %>%
  dplyr::mutate(
    pct_total = 100 * N_total / Biomes_m,
    pct_rank1 = 100 * N_rank1 / Biomes_m,
    Month_lab = factor(month.abb[Month], levels = month.abb),
    Label     = paste0(round(pct_total, 0), "% (", round(pct_rank1, 0), "%)"),
    text_col  = ifelse(pct_total <= 50, "black", "white")  # umbral visual
  )

# # --- ordenar predictores (opción %): por suma de % de rank1 a lo largo del año ---
# order_predictors_rank <- conteo_shap %>%
#   dplyr::group_by(RFE) %>%
#   dplyr::summarise(Total_rank1_pct = sum(pct_rank1, na.rm = TRUE), .groups = "drop") %>%
#   dplyr::arrange(Total_rank1_pct) %>%
#   dplyr::pull(RFE)
# 
# conteo_shap <- conteo_shap %>%
#   dplyr::mutate(RFE = factor(RFE, levels = order_predictors_rank))
# --- ordenar predictores por % de presencia (no por rank1) ---
# --- ordenar predictores por % de presencia (más frecuente arriba) ---
order_predictors_presence <- conteo_shap %>%
  dplyr::group_by(RFE) %>%
  dplyr::summarise(Total_pct_presence = sum(pct_total, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Total_pct_presence)) %>%
  dplyr::pull(RFE)

conteo_shap <- conteo_shap %>%
  dplyr::mutate(RFE = factor(RFE, levels = order_predictors_presence))


# --- heatmap en % ---
p_shap <- ggplot(conteo_shap, aes(x = Month_lab, y = RFE, fill = pct_total)) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = Label, color = text_col), size = 5) +
  scale_color_identity() +
  scale_fill_gradient(limits = c(0, 100), low = "grey90", high = "black", name = "% biomas") +
  labs(
    x = "", y = "Predictors (#1 importance)",
    title = "",
    subtitle = ""
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 14, colour = "black"),
    axis.text.y = element_text(size = 14, colour = "black"),
    axis.title.x = element_text(size = 16, colour = "black", face = "bold"),
    axis.title.y = element_text(size = 16, colour = "black", face = "bold"),
    plot.title   = element_text(size = 24, face = "bold", hjust = 0.5, colour = "black"),
    plot.subtitle= element_text(size = 14, hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = unit(2.5, "cm"),
    legend.title = element_text(size = 14, colour = "black"),
    legend.text  = element_text(size = 14, colour = "black")
  )

print(p_shap)
# Guardado (mismo path que tenías)
ggsave(file.path(output_dir, "matrix_conteo_conteo_predictores_pct.png"),
       plot = p_shap, width = 15, height = 11, dpi = 300)


# ================== LIBRERÍAS ==================
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ================== 1) NORMALIZAR df_shap ==================
# Asume df_shap ya existe con columnas: Biome, Month (opcional), RFE (predictor renombrado), Rank
# Si tu predictor aún se llama 'RFE' tras el recode anterior, seguimos; si no, ajusta el nombre aquí.

# Detectar / crear Year
year_col <- names(df_shap)[grepl("^year$", names(df_shap), ignore.case = TRUE)]
if (length(year_col) == 0) {
  df_shap <- df_shap %>% dplyr::mutate(Year = "")
} else {
  df_shap <- df_shap %>% dplyr::mutate(Year = .data[[year_col[1]]])
  df_shap$Year <- as.character(df_shap$Year)
}

# Rank entero válido
df_shap <- df_shap %>%
  dplyr::mutate(Rank = suppressWarnings(as.integer(Rank))) %>%
  dplyr::filter(!is.na(Rank), Rank >= 1)

# ================== 2) ID de modelo y CONTEOS ==================
has_month <- "Month" %in% names(df_shap)

# Numerador: cuántas veces (modelo único) cada RFE ocupa cada Rank en cada Year
num_rank <- df_shap %>%
  { if (has_month) dplyr::transmute(., Year, RFE, Rank, model_id = paste(Biome, Month, sep = "_"))
    else            dplyr::transmute(., Year, RFE, Rank, model_id = Biome) } %>%
  dplyr::distinct(Year, RFE, Rank, model_id) %>%
  dplyr::group_by(Year, RFE, Rank) %>%
  dplyr::summarise(N_rank = dplyr::n(), .groups = "drop")

# (Opcional) Denominador por si luego quieres %; no se usa en el plot de conteos
denom_rank <- df_shap %>%
  { if (has_month) dplyr::transmute(., Year, Rank, model_id = paste(Biome, Month, sep = "_"))
    else            dplyr::transmute(., Year, Rank, model_id = Biome) } %>%
  dplyr::distinct() %>%
  dplyr::count(Year, Rank, name = "Total_models_rank")

# ================== 3) CUADRÍCULA COMPLETA Y TABLA FINAL ==================
all_ranks <- sort(unique(df_shap$Rank))
all_rfe   <- sort(unique(df_shap$RFE))
all_years <- sort(unique(df_shap$Year))

grid_full <- tidyr::expand_grid(Year = all_years, RFE = all_rfe, Rank = all_ranks)

conteo_year_rank <- grid_full %>%
  dplyr::left_join(num_rank,  by = c("Year", "RFE", "Rank")) %>%
  dplyr::left_join(denom_rank, by = c("Year", "Rank")) %>%
  dplyr::mutate(
    N_rank = tidyr::replace_na(N_rank, 0L),
    Total_models_rank = tidyr::replace_na(Total_models_rank, 0L)
  )

# ================== 4) ORDEN DE PREDICTORES Y ETIQUETAS ==================
# Orden por total de veces que fueron Rank=1 (conteo absoluto)
order_rfe <- conteo_year_rank %>%
  dplyr::filter(Rank == 1) %>%
  dplyr::group_by(RFE) %>%
  dplyr::summarise(TotalRank1 = sum(N_rank, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(TotalRank1)) %>%
  dplyr::pull(RFE)
# 
# conteo_year_rank <- conteo_year_rank %>%
#   dplyr::mutate(
#     RFE    = factor(RFE, levels = order_rfe),
#     Rank_f = factor(Rank, levels = sort(all_ranks)),
#     Label  = ifelse(N_rank > 0, as.character(N_rank), ""),
#     # ----- Color de texto según umbral de CONTEO -----
#     text_col = ifelse(N_rank >= 100, "white", "black")  # <--- UMBRAL EN 100
#   )
# 
# # ================== 5) HEATMAP POR CONTEOS ==================
# p_annual_counts <- ggplot2::ggplot(conteo_year_rank, ggplot2::aes(x = Rank_f, y = RFE, fill = N_rank)) +
#   ggplot2::geom_tile(color = "grey30") +
#   ggplot2::geom_text(ggplot2::aes(label = Label, color = text_col), size = 4) +
#   ggplot2::scale_color_identity() +
#   ggplot2::scale_fill_gradient(low = "grey90", high = "black", name = "# modelos") +
#   ggplot2::labs(
#     x = "Ranking position",
#     y = "",
#     title = "Annual count of occurrences of each predictor at each ranking position"
#     # subtitle = ""
#   ) +
#   ggplot2::facet_wrap(~ Year, scales = "free_y") +
#   ggplot2::theme_minimal() +
#   ggplot2::theme(
#     panel.grid = ggplot2::element_blank(),
#     axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
#     axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
#     axis.title.x = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
#     axis.title.y = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
#     plot.title   = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
#     plot.subtitle= ggplot2::element_text(size = 12, hjust = 0.5, colour = "black"),
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.key.width = grid::unit(2.5, "cm"),
#     legend.title = ggplot2::element_text(size = 12, colour = "black"),
#     legend.text  = ggplot2::element_text(size = 12, colour = "black")
#   ) +
#   ggplot2::scale_y_discrete(drop = FALSE)  # mantiene todos los predictores en todas las facetas
# ... (tu código previo idéntico)

conteo_year_rank <- conteo_year_rank %>%
  dplyr::mutate(
    RFE    = factor(RFE, levels = order_rfe),
    Rank_f = factor(Rank, levels = sort(all_ranks)),
    Label  = ifelse(N_rank > 0, as.character(N_rank), ""),
    text_col = ifelse(N_rank >= 100, "white", "black"),
    N_rank_plot = ifelse(N_rank == 0, NA, N_rank)   # <<< para dejar en blanco los ceros
  )

# ================== 5) HEATMAP POR CONTEOS ==================
p_annual_counts <- ggplot2::ggplot(
  conteo_year_rank,
  ggplot2::aes(x = Rank_f, y = RFE, fill = N_rank_plot)   # <<< usar N_rank_plot
) +
  ggplot2::geom_tile(color = "grey30") +
  ggplot2::geom_text(ggplot2::aes(label = Label, color = text_col), size = 4) +
  ggplot2::scale_color_identity() +
  ggplot2::scale_fill_gradient(
    low = "grey90", high = "black", name = " ",
    na.value = "white"   # <<< celdas con 0 (NA en el mapa) se verán blancas
  ) +
  ggplot2::labs(
    x = "Ranking position",
    y = "",
    title = "Annual count of occurrences of \neach predictor at each ranking position"
  ) +
  ggplot2::facet_wrap(~ Year, scales = "free_y") +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
    axis.text.y = ggplot2::element_text(size = 12, colour = "black"),
    axis.title.x = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
    # axis.title.y = ggplot2::element_text(size = 14, face = "bold", colour = "black"),
    plot.title   = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
    plot.subtitle= ggplot2::element_text(size = 12, hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(2.5, "cm"),
    legend.title = ggplot2::element_text(size = 12, colour = "black"),
    legend.text  = ggplot2::element_text(size = 12, colour = "black")
  ) +
  ggplot2::scale_y_discrete(drop = FALSE)

print(p_annual_counts)

# print(p_annual_counts)
# ... construyes p_annual_counts_base SIN facet:
p_annual_counts_base <- ggplot2::ggplot(
  conteo_year_rank,
  ggplot2::aes(x = Rank_f, y = RFE, fill = N_rank_plot)
) +
  geom_tile(color = "grey30") +
  geom_text(aes(label = Label, color = text_col), size = 4) +
  scale_color_identity() +
  scale_fill_gradient(low = "grey90", high = "black", name = " ", na.value = "white") +
  labs(x = "Ranking position", y = "", title = "Annual count of occurrences of \neach predictor at each ranking position") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5, colour = "black"),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.key.width = grid::unit(2.5, "cm")
  ) +
  scale_y_discrete(drop = FALSE)

# Decidir si facetear:
years_real <- unique(conteo_year_rank$Year[!(conteo_year_rank$Year %in% c("", "All", "ALL"))])

p_annual_counts <- if (length(years_real) > 1) {
  p_annual_counts_base + facet_wrap(~ Year, scales = "free_y")
} else {
  p_annual_counts_base
}

print(p_annual_counts)

# ================== 6) GUARDADO ==================
# Ajusta 'output_dir' a tu ruta; existe en tu flujo anterior
ggplot2::ggsave("/mnt/disco6tb/MRBA60-2/results/D1-Plots/Figure_S3_R1.jpeg",
                plot = p_annual_counts, width = 9, height = 7, dpi = 300)


# ================== 7) PORCENTAJE DE ENTRADA DE CADA PREDICTOR ==================

total_modelos <- 579

tabla_pct_predictores <- conteo_year_rank %>%
  dplyr::group_by(RFE) %>%
  dplyr::summarise(
    N_total = sum(N_rank, na.rm = TRUE),
    Percent_total = round(100 * N_total / 579, 2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(Percent_total))

print(tabla_pct_predictores)



# ==== (Opcional) exportar ====
# write.csv(conteo_shap, file.path(dir_csv, "conteo_shap_summary.csv"), row.names = FALSE)
# ggsave(file.path(dir_csv, "shap_summary_heatmap_rank1.png"), p_shap, width = 10, height = 12, dpi = 300)

# 
# 
# # ==== Paquetes ====
# library(dplyr)
# library(tidyr)
# library(stringr)
# library(ggplot2)
# 
# # ==== Directorio de trabajo ====
# dir_csv <- "/mnt/disco6tb/Dropbox/UAH/FireCCI60/results/M2_RF-2003-2022-CL60-R30/csv/"
# 
# # ==== Diccionario nombres bonitos ====
# rename_map <- c(
#   "f5"               = "FireCCI51",
#   "count_ActiveFire" = "Active Fire",
#   "prec"             = "Precipitation",
#   "temp"             = "Temperature",
#   "FRPsum"           = "FRP sum",
#   "FRPmedian"        = "FRP median",
#   "FWI"              = "FWId95",
#   "lat"              = "Latitude",
#   "LAIHIGH"          = "LAI High",
#   "LAILOW"           = "LAI Low",
#   "wind"             = "Wind Speed",
#   "NDVI"             = "NDVI",
#   "pop"              = "GlobPOP",
#   "cloud"            = "Cloud",
#   "gpp"              = "GPP",
#   "vpd"              = "VPD",
#   "soil"             = "SOIL"
# )
# read_shap_from_csv <- function(fpath) {
#   df <- read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
#   
#   # Arreglar nombres vacíos
#   names(df) <- trimws(names(df))
#   empty_names <- which(is.na(names(df)) | names(df) == "")
#   if (length(empty_names) > 0) {
#     names(df)[empty_names] <- paste0("X", empty_names)
#   }
#   
#   # Detectar columnas (case-insensitive)
#   month_col <- names(df)[grepl("^month$", names(df), ignore.case = TRUE)]
#   var_col   <- names(df)[grepl("variable|feature|predictor", names(df), ignore.case = TRUE)]
#   rank_col  <- names(df)[grepl("rank", names(df), ignore.case = TRUE)]
#   
#   if (length(month_col) == 0 || length(var_col) == 0 || length(rank_col) == 0) {
#     warning("Saltando archivo sin Month/Variable/Rank: ", basename(fpath))
#     return(NULL)
#   }
#   
#   month_col <- month_col[1]; var_col <- var_col[1]; rank_col <- rank_col[1]
#   
#   # Nombre de bioma desde el archivo
#   bio <- basename(fpath)
#   bio <- sub("^shap_summary_(.*)\\.csv$", "\\1", bio)
#   bio <- gsub("__", " ", bio)
#   bio <- gsub("_",  " ", bio)
#   bio <- trimws(bio)
#   
#   df2 <- df %>%
#     dplyr::transmute(
#       Month = suppressWarnings(as.integer(.data[[month_col]])),
#       RFE   = as.character(.data[[var_col]]),
#       Rank  = suppressWarnings(as.integer(.data[[rank_col]])),
#       Biome = bio
#     ) %>%
#     dplyr::filter(!is.na(Month), Month >= 1, Month <= 12)
#   
#   if (nrow(df2) == 0) return(NULL)
#   df2
# }
# 
# 
# # ==== Leer todos los shap_summary ====
# files_shap <- list.files(dir_csv, pattern = "^shap_summary_.*\\.csv$", full.names = TRUE)
# lst_shap <- lapply(files_shap, read_shap_from_csv)
# lst_shap <- lst_shap[!vapply(lst_shap, is.null, logical(1))]
# df_shap <- dplyr::bind_rows(lst_shap)
# 
# # ==== Renombrar predictores ====
# df_shap <- df_shap %>%
#   dplyr::mutate(RFE = dplyr::recode(RFE, !!!rename_map, .default = RFE))
# 
# # ==== Conteo: total y Rank = 1 ====
# conteo_shap <- df_shap %>%
#   dplyr::group_by(Month, RFE) %>%
#   dplyr::summarise(
#     N_total = dplyr::n_distinct(Biome),          # en cuántos biomas entra
#     N_rank1 = sum(Rank == 1, na.rm = TRUE),      # en cuántos fue #1
#     .groups = "drop"
#   ) %>%
#   dplyr::mutate(
#     Month_lab = factor(month.abb[Month], levels = month.abb),
#     Label = paste0(N_total, " (", N_rank1, ")"), # texto combinado
#     text_col = ifelse(N_total <= 25, "black", "white") # color del texto
#   )
# 
# # ==== Ordenar predictores de menor a mayor aportación total ====
# order_predictors <- conteo_shap %>%
#   dplyr::group_by(RFE) %>%
#   dplyr::summarise(Total = sum(N_total), .groups = "drop") %>%
#   dplyr::arrange(Total) %>%
#   dplyr::pull(RFE)
# 
# conteo_shap <- conteo_shap %>%
#   dplyr::mutate(RFE = factor(RFE, levels = order_predictors))
# 
# p_shap <- ggplot(conteo_shap, aes(x = Month_lab, y = RFE, fill = N_total)) +
#   geom_tile(color = "grey30") +
#   geom_text(aes(label = Label, color = text_col), size = 4) +   # ← más pequeño que 5
#   scale_color_identity() +
#   scale_fill_gradient(low = "grey90", high = "black", name = "") +
#   labs(
#     x = "Month", y = "Predictors",
#     title = "Predictor frequency and top-rank importance by month",
#     subtitle = "Total count (times ranked 1)"
#   ) +
#   theme_minimal() +
#   theme(
#     panel.grid = element_blank(),
#     axis.text.x = element_text(size = 14, colour = "black", angle = 0, hjust = 0.5),
#     axis.text.y = element_text(size = 14, colour = "black"),
#     axis.title.x = element_text(size = 16, colour = "black", face = "bold"),
#     axis.title.y = element_text(size = 16, colour = "black", face = "bold"),
#     plot.title = element_text(size = 21, face = "bold", hjust = 0.5, colour = "black"),
#     plot.subtitle = element_text(size = 18, hjust = 0.5, colour = "black"),
#     legend.position = "bottom",
#     legend.direction = "horizontal",
#     legend.key.width = unit(2.5, "cm"),
#     legend.title = element_text(size = 14, colour = "black"),
#     legend.text  = element_text(size = 12, colour = "black")
#     # )
#     # panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
#     # legend.background = element_rect(color = "black", fill = NA, linewidth = 0.8)
#     # legend.box.background = element_rect(color = "black", fill = NA, linewidth = 0.8)
#   )
# 
# print(p_shap)
# 
# # ==== (Opcional) exportar ====
# # write.csv(conteo_wide, file.path(dir_csv, "predictor_month_counts_wide.csv"), row.names = FALSE)
# output_dir= "/mnt/disco6tb/Dropbox/UAH/FireCCI60/ATBD"
# ggsave(file.path(output_dir, "matrix_conteo_conteo_predictores.png"),
#        plot = p_shap, width = 10, height = 12, dpi = 300)

# ==== (Opcional) exportar ====
# write.csv(conteo_shap, file.path(dir_csv, "conteo_shap_summary.csv"), row.names = FALSE)
# ggsave(file.path(dir_csv, "shap_summary_heatmap.png"), p_shap, width = 10, height = 12, dpi = 300)

