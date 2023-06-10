# calculate rmse
rmse <- function(y, y_pred) {
  (mean((y - y_pred)^2))^0.5
}
