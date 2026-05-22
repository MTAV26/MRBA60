build_prediction_df <- function(arr, tt, idx, include_response = FALSE) {
  
  out <- data.frame(
    f5 = arr$f51[, , tt][idx],
    count_ActiveFire = arr$af[, , tt][idx],
    prec = arr$prec[, , tt][idx],
    temp = arr$temp[, , tt][idx],
    FRPsum = arr$frp_sum[, , tt][idx],
    FRPmedian = arr$frp_median[, , tt][idx],
    NDVI = arr$ndvi[, , tt][idx],
    FWI = arr$fwi[, , tt][idx],
    wind = arr$wind[, , tt][idx],
    lat = arr$lat[, , tt][idx],
    lon = arr$lon[, , tt][idx],
    cloud = arr$cloud[, , tt][idx],
    vpd = arr$vpd[, , tt][idx],
    soil = arr$soil[, , tt][idx]
  )
  
  if (include_response) {
    out <- cbind(
      f3 = arr$s3[, , tt][idx],
      out
    )
  }
  
  out
}
