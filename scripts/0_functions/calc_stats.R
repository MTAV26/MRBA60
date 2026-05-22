calc_stats <- function(y_true, y_pred) {
  
  ok <- is.finite(y_true) & is.finite(y_pred)
  
  if (sum(ok) < 2) {
    return(data.frame(
      bias = NA_real_,
      NME = NA_real_,
      RMSE = NA_real_,
      R2 = NA_real_,
      Correlacion = NA_real_,
      p_value = NA_real_,
      Significancia = NA_character_
    ))
  }
  
  y_true <- y_true[ok]
  y_pred <- y_pred[ok]
  
  denom_nme <- sum(abs(y_true - mean(y_true)))
  denom_r2 <- sum((y_true - mean(y_true))^2)
  
  bias <- mean(y_pred - y_true)
  NME <- ifelse(denom_nme > 0, sum(abs(y_true - y_pred)) / denom_nme, NA_real_)
  RMSE <- sqrt(mean((y_pred - y_true)^2))
  R2 <- ifelse(denom_r2 > 0, 1 - sum((y_true - y_pred)^2) / denom_r2, NA_real_)
  
  cor_test <- tryCatch(
    suppressWarnings(stats::cor.test(y_true, y_pred, method = "spearman")),
    error = function(e) NULL
  )
  
  if (is.null(cor_test)) {
    Correlacion <- NA_real_
    p_value <- NA_real_
  } else {
    Correlacion <- as.numeric(cor_test$estimate)
    p_value <- cor_test$p.value
  }
  
  Significancia <- ifelse(is.finite(p_value) && p_value < 0.05, "*", "ns")
  
  data.frame(
    bias = bias,
    NME = NME,
    RMSE = RMSE,
    R2 = R2,
    Correlacion = Correlacion,
    p_value = p_value,
    Significancia = Significancia
  )
}