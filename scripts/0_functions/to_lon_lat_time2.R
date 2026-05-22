
to_lon_lat_time2 <- function(x, lon, lat) {
  
  dx <- dim(x)
  
  if (length(dx) != 3) {
    stop("El objeto no tiene 3 dimensiones.")
  }
  
  if (dx[1] == length(lon) && dx[2] == length(lat)) {
    return(x)
  }
  
  if (dx[1] == length(lat) && dx[2] == length(lon)) {
    return(aperm(x, c(2, 1, 3)))
  }
  
  stop(
    "Las dimensiones espaciales no coinciden con lon/lat.\n",
    "dim = ", paste(dx, collapse = " x "), "\n",
    "length(lon) = ", length(lon), "\n",
    "length(lat) = ", length(lat)
  )
}
