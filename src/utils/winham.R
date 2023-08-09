source("src/utils/mse.R")
winham <- function(train_data,
                   test_data,
                   rf_model,
                   norm_param,
                   rf_predictions_train_all,
                   rf_predictions_test_all) {
  inbcs <- do.call(cbind, rf_model$inbag.counts)
  rf_predictions_train_all[which(inbcs > 0)] <- NA
  weights_unnormalized <- colMeans(
    abs(rf_predictions_train_all - train_data$target),
    na.rm = TRUE
  )
  weights <- (exp(1 / weights_unnormalized) / sum(exp(1 / weights_unnormalized)))
  preds_winham <- rf_predictions_test_all %*% weights
  preds_winham <- (preds_winham * norm_param$sd) + norm_param$mean
  mse_winham <- mse(preds_winham, test_data$target)
  return(mse_winham)
}
