build_monthly_training_data <- function(arr_common, mes) {
  
  mes_indices_central <- which(lubridate::month(dates_common) == mes)
  
  mes_indices_window <- sort(unique(c(
    mes_indices_central,
    mes_indices_central - 1,
    mes_indices_central + 1
  )))
  
  mes_indices_window <- mes_indices_window[
    mes_indices_window >= 1 &
      mes_indices_window <= length(dates_common)
  ]
  
  monthly_data <- list()
  
  for (tt in mes_indices_window) {
    
    f3_layer <- arr_common$s3[, , tt]
    f5_layer <- arr_common$f51[, , tt]
    af_layer <- arr_common$af[, , tt]
    
    idx <- which(
      is.finite(f3_layer) &
        is.finite(f5_layer) &
        is.finite(af_layer)
    )
    
    if (length(idx) > 0) {
      
      current_year <- lubridate::year(dates_common[tt])
      key <- as.character(current_year)
      
      df_new <- build_prediction_df(
        arr = arr_common,
        tt = tt,
        idx = idx,
        include_response = TRUE
      )
      
      df_new$year <- current_year
      df_new$month <- lubridate::month(dates_common[tt])
      
      if (!is.null(monthly_data[[key]])) {
        monthly_data[[key]] <- rbind(monthly_data[[key]], df_new)
      } else {
        monthly_data[[key]] <- df_new
      }
    }
  }
  
  if (length(monthly_data) == 0) {
    return(list(
      df_full = NULL,
      mes_indices_central = mes_indices_central,
      mes_indices_window = mes_indices_window
    ))
  }
  
  df_full <- do.call(rbind, monthly_data)
  
  df_full <- df_full[
    complete.cases(df_full[, c(response, candidate_predictors), drop = FALSE]),
  ]
  
  list(
    df_full = df_full,
    mes_indices_central = mes_indices_central,
    mes_indices_window = mes_indices_window
  )
}