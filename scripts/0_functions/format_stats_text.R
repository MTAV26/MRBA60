format_stats_text <- function(stats_tbl) {
  
  if (is.null(stats_tbl) || nrow(stats_tbl) == 0) {
    return(NULL)
  }
  
  lines <- "Performance Comparisons:"
  
  for (ii in seq_len(nrow(stats_tbl))) {
    model_name <- as.character(stats_tbl$Modelo[ii])
    sig <- ""
    
    if ("Significancia" %in% names(stats_tbl)) {
      sig <- as.character(stats_tbl$Significancia[ii])
      if (is.na(sig) || !nzchar(sig)) {
        sig <- ""
      }
    }
    
    lines <- c(
      lines,
      model_name,
      paste0(
        "Bias = ", format_stat_value(stats_tbl$bias[ii]),
        " | NME = ", format_stat_value(stats_tbl$NME[ii]),
        " | RMSE = ", format_stat_value(stats_tbl$RMSE[ii])
      ),
      paste0(
        "R2 = ", format_stat_value(stats_tbl$R2[ii]),
        " | Cor = ", format_stat_value(stats_tbl$Correlacion[ii]), sig
      )
    )
    
    if (ii < nrow(stats_tbl)) {
      lines <- c(lines, "")
    }
  }
  
  paste(lines, collapse = "\n")
}