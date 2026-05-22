predict_layer_rf <- function(modelo, arr, tt, best_preds, area_matrix_crop, base_layer) {
  
  layer <- base_layer
  
  f5_layer <- arr$f51[, , tt]
  af_layer <- arr$af[, , tt]
  
  idx <- which(
    is.finite(f5_layer) &
      is.finite(af_layer)
  )
  
  if (length(idx) == 0) {
    return(layer)
  }
  
  pred_df <- build_prediction_df(
    arr = arr,
    tt = tt,
    idx = idx,
    include_response = FALSE
  )
  
  ok <- complete.cases(pred_df[, best_preds, drop = FALSE])
  
  if (!any(ok)) {
    return(layer)
  }
  
  pred_df <- pred_df[ok, , drop = FALSE]
  idx_ok <- idx[ok]
  
  preds <- predict(modelo, newdata = pred_df)
  cap_phys <- area_matrix_crop[idx_ok]
  
  preds <- pmax(preds, 0)
  preds <- pmin(preds, cap_phys)
  
  layer[idx_ok] <- preds
  layer
}