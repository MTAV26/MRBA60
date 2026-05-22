safe_name <- function(x) {
  x <- gsub(" ", "_", x)
  x <- gsub("[^[:alnum:]_]", "", x)
  x
}
