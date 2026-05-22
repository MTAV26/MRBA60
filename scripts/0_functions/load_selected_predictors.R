load_selected_predictors <- function(safe_biome_name) {
  
  csv_file <- file.path(
    selected_predictors_dir2,
    paste0("SelectedPredictors_", safe_biome_name, "_COMMON.csv")
  )
  
  if (!file.exists(csv_file)) {
    return(NULL)
  }
  
  tbl <- read.csv(csv_file, stringsAsFactors = FALSE)
  
  if (!"Month" %in% names(tbl)) {
    stop("El CSV de predictores no contiene la columna Month: ", csv_file)
  }
  
  predictor_col <- NA_character_
  formula_col <- NA_character_
  
  if ("Predictors_after_RFE" %in% names(tbl)) {
    predictor_col <- "Predictors_after_RFE"
  } else if ("RFE_selected" %in% names(tbl)) {
    predictor_col <- "RFE_selected"
  }
  
  if ("Formula_after_RFE" %in% names(tbl)) {
    formula_col <- "Formula_after_RFE"
  } else if ("Formula" %in% names(tbl)) {
    formula_col <- "Formula"
  }
  
  by_month <- vector("list", 12)
  
  for (mes in 1:12) {
    row_mes <- tbl[tbl$Month == mes, , drop = FALSE]
    
    if (nrow(row_mes) == 0) {
      next
    }
    
    row_mes <- row_mes[nrow(row_mes), , drop = FALSE]
    preds <- character(0)
    
    if (!is.na(predictor_col)) {
      preds <- parse_predictor_string(row_mes[[predictor_col]])
    }
    
    if (length(preds) == 0 && !is.na(formula_col)) {
      preds <- parse_predictor_string(row_mes[[formula_col]])
    }
    
    by_month[[mes]] <- preds
  }
  
  list(
    file = csv_file,
    table = tbl,
    by_month = by_month
  )
}