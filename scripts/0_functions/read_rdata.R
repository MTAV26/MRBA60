read_rdata <- function(path, object = NULL) {
  e <- new.env(parent = emptyenv())
  nm <- load(path, envir = e)
  
  if (!is.null(object)) {
    if (!object %in% nm) {
      stop(
        "El objeto '", object, "' no existe en ", basename(path),
        ". Objetos disponibles: ", paste(nm, collapse = ", ")
      )
    }
    return(e[[object]])
  }
  
  if (length(nm) == 1) {
    return(e[[nm]])
  }
  
  stop(
    "El archivo ", basename(path), " contiene varios objetos: ",
    paste(nm, collapse = ", "),
    ". Indica explícitamente cuál cargar con object = ..."
  )
}
