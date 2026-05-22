read_nc_var <- function(fname, varname) {
  nc <- ncdf4::nc_open(fname)
  on.exit(ncdf4::nc_close(nc), add = TRUE)
  ncdf4::ncvar_get(nc, varname)
}