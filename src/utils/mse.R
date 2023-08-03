# calculate rmse
mse <- function(y, y_pred) {
  mean((y - y_pred)^2)
}
