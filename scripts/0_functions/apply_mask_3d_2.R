apply_mask_3d_2 <- function(arr, mask_lonlat) {
  stopifnot(length(dim(arr)) == 3)
  stopifnot(all(dim(arr)[1:2] == dim(mask_lonlat)))
  
  for (tt in seq_len(dim(arr)[3])) {
    z <- arr[, , tt]
    z[!mask_lonlat] <- NA
    arr[, , tt] <- z
  }
  
  arr
}