
format_stat_value <- function(x, digits = 2) {
  x <- suppressWarnings(as.numeric(x[1]))
  
  if (!is.finite(x)) {
    return("NA")
  }
  
  if (abs(x) >= 100) {
    return(formatC(x, format = "f", digits = 0, big.mark = ","))
  }
  
  formatC(x, format = "f", digits = digits)
}
