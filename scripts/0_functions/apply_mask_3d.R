apply_mask_3d <- function(arr, mask_matrix) {
  out <- arr
  for (tt in seq_len(dim(out)[3])) {
    layer <- out[, , tt]
    layer[!mask_matrix] <- NA
    out[, , tt] <- layer
  }
  out
}