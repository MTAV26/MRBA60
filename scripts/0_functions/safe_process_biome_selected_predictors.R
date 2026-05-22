
safe_process_biome_selected_predictors <- function(bioma) {
  
  safe_biome_name <- safe_name(as.character(bioma))
  
  tryCatch(
    {
      process_biome_selected_predictors(bioma)
    },
    error = function(e) {
      
      error_file <- file.path(
        output_dir_log,
        paste0("ERROR_SelectedPredictors_", safe_biome_name, ".txt")
      )
      
      msg <- paste0(
        "ERROR EN BIOMA: ", bioma, "\n",
        "Fecha: ", as.character(Sys.time()), "\n",
        "Mensaje: ", conditionMessage(e), "\n"
      )
      
      writeLines(msg, error_file)
      
      list(
        Biome = safe_biome_name,
        Status = "ERROR",
        n_months_selected = NA_integer_,
        seconds = NA_real_,
        csv = NA_character_,
        error = conditionMessage(e)
      )
    }
  )
}
