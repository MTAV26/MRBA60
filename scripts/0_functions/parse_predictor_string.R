parse_predictor_string <- function(x) {
  
  if (length(x) == 0 || is.na(x) || !nzchar(trimws(x))) {
    return(character(0))
  }
  
  x <- as.character(x[1])
  x <- gsub("^.*~", "", x)
  x <- gsub("`", "", x)
  
  out <- unlist(strsplit(x, "\\+"))
  out <- trimws(out)
  out <- out[nzchar(out)]
  out <- unique(out)
  out <- out[!tolower(out) %in% "gpp"]
  out
}