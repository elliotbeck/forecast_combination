source("src/utils/get_mse.R")
cesaro <- function(train_data,
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
  weights_position <- sort(weights_unnormalized, index.return = TRUE, decreasing = TRUE)$ix
  weights <- cumsum(rep(1 / length(weights_unnormalized), length(weights_unnormalized)))
  weights <- cbind(weights, weights_position)
  weights <- weights[order(weights[, 2]), 1]
  weights <- weights / sum(weights)
  preds_cesaro <- rf_predictions_test_all %*% weights
  preds_cesaro <- (preds_cesaro * norm_param$sd) + norm_param$mean
  mse_cesaro <- mse(preds_cesaro, test_data$target)
  return(mse_cesaro)
}
