rm(list = ls())
graphics.off()
gc()

# =========================================================
# Librerías
# =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(tibble)
})

# =========================================================
# 1) Directorio de los CSV
# =========================================================

dir_csv <- "/mnt/disco6tb/MRBA60-2/results/B1-MRBA60-2003-2024/csv"

# =========================================================
# 2) Listar archivos SelectedPredictors_*.csv
# =========================================================

files <- list.files(
  dir_csv,
  pattern = "^SelectedPredictors_.*\\.csv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No SelectedPredictors_*.csv files found. Check dir_csv.")
}

cat("Number of SelectedPredictors files found:", length(files), "\n")

# =========================================================
# 3) Función para limpiar nombres de biomas
# =========================================================

clean_biome_name <- function(x) {
  
  bio <- x
  
  # Eliminar extensión y prefijo del archivo
  bio <- basename(bio)
  bio <- sub("^SelectedPredictors_(.*)\\.csv$", "\\1", bio)
  
  # ---------------------------------------------------------
  # Limpieza de tokens que no pertenecen al nombre del bioma
  # ---------------------------------------------------------
  
  # Eliminar Used si aparece como token independiente:
  # Used_Biome, Biome_Used, Biome__Used__, etc.
  bio <- gsub("(^|_+)USED(_+|$)", "_", bio, ignore.case = TRUE)
  
  # Eliminar COMMON si aparece como token independiente:
  # COMMON_Biome, Biome_COMMON, Biome__COMMON__, etc.
  bio <- gsub("(^|_+)COMMON(_+|$)", "_", bio, ignore.case = TRUE)
  
  # Limpieza adicional por seguridad
  bio <- gsub("^USED_+", "", bio, ignore.case = TRUE)
  bio <- gsub("_+USED$", "", bio, ignore.case = TRUE)
  bio <- gsub("^COMMON_+", "", bio, ignore.case = TRUE)
  bio <- gsub("_+COMMON$", "", bio, ignore.case = TRUE)
  
  # Si quedara algún token aislado tras limpiezas anteriores
  bio <- gsub("\\bUSED\\b", "", bio, ignore.case = TRUE)
  bio <- gsub("\\bCOMMON\\b", "", bio, ignore.case = TRUE)
  
  # Convertir guiones bajos en espacios
  bio <- gsub("_+", " ", bio)
  
  # Limpiar espacios duplicados, iniciales y finales
  bio <- gsub("\\s+", " ", bio)
  bio <- trimws(bio)
  
  return(bio)
}

# =========================================================
# 4) Leer columna Month y extraer nombre limpio del bioma
# =========================================================

read_months_from_csv <- function(fpath) {
  
  df <- read.csv(
    fpath,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  # Detectar columna Month ignorando mayúsculas/minúsculas
  month_col <- names(df)[grepl("^month$", names(df), ignore.case = TRUE)]
  
  if (length(month_col) == 0) {
    warning("No Month column found in: ", basename(fpath))
    return(NULL)
  }
  
  # Extraer meses válidos
  months_ok <- unique(df[[month_col[1]]])
  months_ok <- suppressWarnings(as.integer(months_ok))
  months_ok <- months_ok[!is.na(months_ok) & months_ok >= 1 & months_ok <= 12]
  
  if (length(months_ok) == 0) {
    warning("No valid months found in: ", basename(fpath))
    return(NULL)
  }
  
  # Extraer bioma limpio desde el nombre del archivo
  bio <- clean_biome_name(fpath)
  
  tibble(
    Biome = bio,
    Month = months_ok,
    Calculated = TRUE
  )
}

# =========================================================
# 5) Aplicar lectura a todos los CSV
# =========================================================

lst <- lapply(files, read_months_from_csv)
lst <- lst[!vapply(lst, is.null, logical(1))]

if (length(lst) == 0) {
  stop("No valid SelectedPredictors files with a Month column were found.")
}

df_months <- bind_rows(lst)

# =========================================================
# 6) Comprobaciones
# =========================================================

cat("Number of biomes detected:", length(unique(df_months$Biome)), "\n")
cat("Biome names:\n")
print(sort(unique(df_months$Biome)))

# Comprobación explícita de que Used no queda en el eje Y
if (any(grepl("\\bUsed\\b", df_months$Biome, ignore.case = TRUE))) {
  warning("Some biome names still contain 'Used'. Check file naming pattern.")
  print(sort(unique(df_months$Biome[grepl("\\bUsed\\b", df_months$Biome, ignore.case = TRUE)])))
}

# =========================================================
# 7) Completar meses faltantes 1..12 y ordenar biomas
# =========================================================

df_summary <- df_months %>%
  group_by(Biome) %>%
  complete(Month = 1:12) %>%
  mutate(
    Calculated = ifelse(is.na(Calculated), FALSE, Calculated)
  ) %>%
  ungroup() %>%
  mutate(
    Month_lab = factor(month.abb[Month], levels = month.abb),
    Biome = factor(Biome, levels = sort(unique(Biome)))
  )

# =========================================================
# 8) Plot heatmap Bioma/Mes
# =========================================================

p_calc <- ggplot(df_summary, aes(x = Month_lab, y = fct_rev(Biome))) +
  geom_tile(
    aes(fill = Calculated),
    color = "black",
    linewidth = 0.25
  ) +
  scale_fill_manual(
    values = c(
      "FALSE" = "#DEAB76",
      "TRUE"  = "#76DEC2"
    ),
    labels = c(
      "FALSE" = "Not calculated",
      "TRUE"  = "Calculated"
    ),
    drop = FALSE
  ) +
  labs(
    x = "Month",
    y = "Biomes",
    title = "Harmonised Biomes/Month",
    subtitle = "(≥ 30 valid grid cells in ≥ 2 years)",
    fill = "Harmonised"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5)
  )

print(p_calc)

# =========================================================
# 9) Guardar figura
# =========================================================

output_dir <- "/mnt/disco6tb/MRBA60-2/results/D1-Plots/"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

ggsave(
  filename = file.path(output_dir, "FigureS2_R1.png"),
  plot = p_calc,
  width = 10,
  height = 12,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "FigureS2_R1.pdf"),
  plot = p_calc,
  width = 10,
  height = 12
)

cat("Figures saved in:\n")
cat(file.path(output_dir, "FigureS2_R1.png"), "\n")
cat(file.path(output_dir, "FigureS2_R1.pdf"), "\n")
