
to_lon_lat_time <- function(x, grid_order) {
  if (length(dim(x)) != 3) {
    stop("El objeto no tiene 3 dimensiones.")
  }
  
  if (grid_order == "lat_lon_time") {
    return(aperm(x, c(2, 1, 3)))
  } else if (grid_order == "lon_lat_time") {
    return(x)
  } else {
    stop("grid_order no reconocido: ", grid_order)
  }
}