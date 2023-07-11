# Libraries
library(CVXR)
source("src/utils/rmse.R")

get_performance <- function(
    mean_vector,
    cov_matrix,
    predictions_train,
    predictions_test,
    labels_train,
    labels_test,
    kappa) {
  w <- Variable(ncol(cov_matrix))
  objective <- square((t(w) %*% mean_vector)) + quad_form(w, cov_matrix)
  constraints <- list(
    sum(w) == 1,
    sum(abs(w)) <= kappa
  )
  prob <- Problem(Minimize(objective), constraints)
  solution <- solve(prob)

  # Weighted prediction
  predictions_weighted_train <- as.matrix(predictions_train) %*% solution$getValue(w)
  predictions_weighted_test <- as.matrix(predictions_test) %*% solution$getValue(w)

  # Calculate rmse
  rmse_train <- rmse(predictions_weighted_train, labels_train)
  rmse_test <- rmse(predictions_weighted_test, labels_test)

  return(c(rmse_train, rmse_test))
}
