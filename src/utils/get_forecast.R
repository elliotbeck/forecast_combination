# Libraries
source("src/utils/get_weights.R")

get_forecast <- function(
    mean_vector,
    cov_matrix,
    predictions_test,
    kappa = 2,
    mean = 0,
    sd = 1) {
    # Get weights
    w <- get_weights(mean_vector, cov_matrix, kappa)

    # Weighted prediction
    predictions_weighted_test <- as.matrix(predictions_test) %*% w
    predictions_weighted_test <- (predictions_weighted_test * sd) + mean
    return(predictions_weighted_test)
}
