rm(list = ls())
graphics.off()
gc()

# =========================================================
# LIBRERÍAS
# =========================================================

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(tibble)
  library(purrr)
})

# =========================================================
# 1) DIRECTORIO Y ARCHIVOS
# =========================================================

dir_csv <- "/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/csv"

files <- list.files(
  dir_csv,
  pattern = "^SelectedPredictors_.*\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No se encontraron archivos SelectedPredictors_*.csv en dir_csv.")
}

cat("Número de archivos encontrados:", length(files), "\n")

# =========================================================
# 2) CONFIGURACIÓN GENERAL
# =========================================================

# Denominador fijo: número total de combinaciones bioma-mes consideradas
N_total <- 579

# Diccionario de nombres bonitos para el plot
rename_map <- c(
  "f5"               = "FireCCI51",
  "count_ActiveFire" = "AF30F",
  "prec"             = "Precipitation mean",
  "temp"             = "Temperature mean",
  "FRPsum"           = "FRP sum",
  "FRPmedian"        = "FRP median",
  "FWI"              = "FWI95d",
  "lat"              = "Latitude",
  "lon"              = "Longitude",
  "wind"             = "Wind Speed",
  "NDVI"             = "NDVI",
  "cloud"            = "Cloud",
  "vpd"              = "VPD",
  "soil"             = "Soil moisture"
)

nice_name <- function(x) {
  dplyr::recode(x, !!!rename_map, .default = x)
}

# Colores
fill_colors <- c(
  "After removing the autocorrelation" = "#DEAB76",
  "After the RFE selected"             = "#76DEC2"
)

# =========================================================
# 3) FUNCIÓN PARA EXTRAER PREDICTORES
# =========================================================
# Sirve tanto para fórmulas:
#   burned_area ~ f5 + NDVI + temp
# como para listas:
#   f5, NDVI, temp
#   f5 + NDVI + temp
#   c("f5", "NDVI", "temp")
# =========================================================

extract_predictors <- function(x) {
  
  if (is.null(x) || all(is.na(x))) return(character(0))
  
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  
  if (length(x) == 0) return(character(0))
  
  # Quitar parte izquierda de fórmula si existe
  x <- sub("^[^~]*~", "", x)
  
  # Quitar sintaxis tipo c(...)
  x <- gsub("c\\(", "", x)
  x <- gsub("\\)", "", x)
  
  # Quitar comillas
  x <- gsub("\"", "", x)
  x <- gsub("'", "", x)
  
  # Separar por +, coma, punto y coma
  tokens <- unlist(strsplit(x, "[+;,]+"))
  tokens <- trimws(tokens)
  tokens <- tokens[nzchar(tokens)]
  
  if (length(tokens) == 0) return(character(0))
  
  # Limpiar funciones tipo log(x), I(x^2), s(x), etc.
  tokens <- gsub("\\bI\\s*\\(([^)]*)\\)", "\\1", tokens)
  tokens <- gsub("\\b(log|sqrt|scale|as.numeric|as.factor|factor)\\s*\\(([^)]*)\\)", "\\2", tokens)
  
  # Quitar espacios
  tokens <- gsub("[[:space:]]", "", tokens)
  
  # Separar interacciones u operadores
  tokens <- gsub("[:*^/|-]", " ", tokens)
  tokens <- trimws(tokens)
  
  # Re-separar por espacios tras limpiar operadores
  tokens <- unlist(strsplit(tokens, "\\s+"))
  tokens <- tokens[nzchar(tokens)]
  
  # Eliminar constante 1
  tokens <- tokens[tokens != "1"]
  
  # Mantener solo nombres válidos
  tokens <- gsub("[^A-Za-z0-9_]", "", tokens)
  tokens <- tokens[nzchar(tokens)]
  
  unique(tokens)
}

# =========================================================
# 4) FUNCIÓN PARA CONTAR PREDICTORES POR ARCHIVO
# =========================================================

count_from_file <- function(f) {
  
  df <- suppressWarnings(
    readr::read_csv(f, show_col_types = FALSE, guess_max = 10000)
  )
  
  if (nrow(df) == 0) return(tibble())
  
  cn <- names(df)
  
  # ---------------------------------------------------------
  # Detectar columna después de eliminar autocorrelación
  # ---------------------------------------------------------
  # Prioridad:
  # 1) columna explícita de predictores
  # 2) columna de fórmula
  # 3) búsqueda flexible por autocorrelation/autocorr/formula
  # ---------------------------------------------------------
  
  col_auto <- cn[grepl("^Predictors_after_autocorrelation$", cn, ignore.case = TRUE)]
  
  if (length(col_auto) == 0) {
    col_auto <- cn[grepl("^Formula_after_autocorrelation$", cn, ignore.case = TRUE)]
  }
  
  if (length(col_auto) == 0) {
    col_auto <- cn[
      grepl("autocorr|autocorrelation", cn, ignore.case = TRUE) &
        !grepl("rfe", cn, ignore.case = TRUE)
    ]
  }
  
  # ---------------------------------------------------------
  # Detectar columna después del RFE
  # ---------------------------------------------------------
  
  col_rfe <- cn[grepl("^Predictors_after_RFE$", cn, ignore.case = TRUE)]
  
  if (length(col_rfe) == 0) {
    col_rfe <- cn[grepl("^Formula_after_RFE$", cn, ignore.case = TRUE)]
  }
  
  if (length(col_rfe) == 0) {
    col_rfe <- cn[grepl("rfe", cn, ignore.case = TRUE)]
  }
  
  if (length(col_auto) == 0 && length(col_rfe) == 0) {
    warning("No se encontraron columnas de predictores en: ", basename(f))
    return(tibble())
  }
  
  # ---------------------------------------------------------
  # Conteo después de eliminar autocorrelación
  # ---------------------------------------------------------
  
  occ_auto <- tibble()
  
  if (length(col_auto) > 0) {
    
    preds_by_row <- lapply(df[[col_auto[1]]], extract_predictors)
    all_preds <- unique(unlist(preds_by_row))
    
    if (length(all_preds) > 0) {
      occ_auto <- tibble(
        File      = basename(f),
        Predictor = all_preds,
        n_present = vapply(
          all_preds,
          function(p) {
            sum(vapply(preds_by_row, function(v) p %in% v, logical(1)))
          },
          integer(1)
        ),
        N = length(preds_by_row),
        Stage = "After removing the autocorrelation"
      )
    }
  }
  
  # ---------------------------------------------------------
  # Conteo después del RFE
  # ---------------------------------------------------------
  
  occ_rfe <- tibble()
  
  if (length(col_rfe) > 0) {
    
    preds_by_row <- lapply(df[[col_rfe[1]]], extract_predictors)
    all_preds <- unique(unlist(preds_by_row))
    
    if (length(all_preds) > 0) {
      occ_rfe <- tibble(
        File      = basename(f),
        Predictor = all_preds,
        n_present = vapply(
          all_preds,
          function(p) {
            sum(vapply(preds_by_row, function(v) p %in% v, logical(1)))
          },
          integer(1)
        ),
        N = length(preds_by_row),
        Stage = "After the RFE selected"
      )
    }
  }
  
  dplyr::bind_rows(occ_auto, occ_rfe)
}

# =========================================================
# 5) LEER TODOS LOS ARCHIVOS Y COMBINAR
# =========================================================

counts_list <- purrr::map(files, count_from_file)
counts <- dplyr::bind_rows(counts_list)

if (nrow(counts) == 0) {
  stop("No se encontraron predictores en las columnas de autocorrelación o RFE.")
}

cat("\nResumen de fases encontradas:\n")
print(table(counts$Stage))

cat("\nPrimeras filas de counts:\n")
print(head(counts))

# =========================================================
# 6) AGREGAR SOBRE TODOS LOS BIOMAS / MESES
# =========================================================

agg <- counts %>%
  dplyr::group_by(Stage, Predictor) %>%
  dplyr::summarise(
    n_present = sum(n_present, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Percent = 100 * n_present / N_total
  )

cat("\nTabla agregada:\n")
print(agg)

# =========================================================
# 7) CREAR TABLA WIDE AUTO / RFE
# =========================================================

wide <- agg %>%
  dplyr::select(Stage, Predictor, n_present) %>%
  tidyr::pivot_wider(
    names_from  = Stage,
    values_from = n_present,
    values_fill = 0
  )

# Crear columnas si falta alguna fase
if (!"After removing the autocorrelation" %in% names(wide)) {
  wide[["After removing the autocorrelation"]] <- 0
}

if (!"After the RFE selected" %in% names(wide)) {
  wide[["After the RFE selected"]] <- 0
}

wide <- wide %>%
  dplyr::rename(
    Auto = `After removing the autocorrelation`,
    RFE  = `After the RFE selected`
  ) %>%
  dplyr::mutate(
    pct_auto = 100 * Auto / N_total,
    pct_rfe  = 100 * RFE  / N_total
  )

# =========================================================
# 8) ORDEN DEL EJE X
# =========================================================
# Se ordena por el porcentaje tras el RFE

order_ref <- wide %>%
  dplyr::arrange(dplyr::desc(pct_rfe)) %>%
  dplyr::pull(Predictor)

order_ref_nice <- nice_name(order_ref)

wide_plot <- wide %>%
  dplyr::mutate(
    Predictor_nice = factor(
      nice_name(Predictor),
      levels = order_ref_nice
    )
  )

# =========================================================
# 9) PLOT CON BARRAS SOLAPADAS (MISMA ANCHURA)
# =========================================================

p_overlap <- ggplot(wide_plot, aes(x = Predictor_nice)) +
  
  geom_col(
    aes(
      y = pct_auto,
      fill = "After removing the autocorrelation"
    ),
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  
  geom_col(
    aes(
      y = pct_rfe,
      fill = "After the RFE selected"
    ),
    width = 0.75,
    color = "black",
    linewidth = 0.25
  ) +
  
  scale_fill_manual(
    values = c(
      "After removing the autocorrelation" = "#DEAB76",
      "After the RFE selected"             = "#76DEC2"
    ),
    breaks = c(
      "After removing the autocorrelation",
      "After the RFE selected"
    ),
    drop = FALSE
  ) +
  
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 10),
    expand = expansion(mult = c(0, 0.01))
  ) +
  
  labs(
    title = paste0("Predictor occurrence % (All biomes/Month; N=", N_total, ")"),
    x = NULL,
    y = "%",
    fill = NULL
  ) +
  
  theme_classic(base_size = 18) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
    axis.title.y = element_text(face = "bold", size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
    axis.text.y = element_text(size = 16),
    legend.position = "top",
    legend.text = element_text(size = 17),
    legend.key.height = unit(16, "pt"),
    legend.key.width  = unit(24, "pt"),
    plot.margin = ggplot2::margin(t = 8, r = 14, b = 28, l = 14)
  )

print(p_overlap)

# =========================================================
# 10) GUARDAR FIGURA
# =========================================================

out_dir <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(out_dir, "Figure_S3_R1_PredictorOccurrence_overlap_Auto_RFE.pdf"),
  plot     = p_overlap,
  width    = 15,
  height   = 8,
  units    = "in"
)

ggsave(
  filename = file.path(out_dir, "Figure_S3_R1_PredictorOccurrence_overlap_Auto_RFE.jpeg"),
  plot     = p_overlap,
  width    = 15,
  height   = 8,
  units    = "in",
  dpi      = 300,
  device   = "jpeg"
)

cat("\nFigura guardada en:\n")
cat(file.path(out_dir, "Figure_S3_R1_PredictorOccurrence_overlap_Auto_RFE.pdf"), "\n")
cat(file.path(out_dir, "Figure_S3_R1_PredictorOccurrence_overlap_Auto_RFE.jpeg"), "\n")

