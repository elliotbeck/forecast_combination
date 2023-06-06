# calculate rmse
rmse <- function(x, cpi) {
  (mean((x - cpi)^2))^0.5
}
