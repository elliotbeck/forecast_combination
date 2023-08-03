# Libraries
source("src/utils/mse.R")
source("src/simulations/get_weights.R")

get_performance <- function(
    mean_vector,
    cov_matrix,
    predictions_test,
    labels_test,
    kappa,
    mean,
    sd) {
  # Get weights
  w <- get_weights(mean_vector, cov_matrix, kappa)

  # Weighted prediction
  predictions_weighted_test <- as.matrix(predictions_test) %*% w
  predictions_weighted_test <- (predictions_weighted_test * sd) + mean

  # Calculate rmse
  rmse_test <- rmse(predictions_weighted_test, labels_test)

  return(rmse_test)
}
