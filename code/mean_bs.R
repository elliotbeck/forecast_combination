mean_bs <- function (x) {
  p <- nrow(x)
  n <- ncol(x)
  cov_mtrx <- qis(t(x)) 
  invSS <- solve(cov_mtrx)
  means <- .rowMeans(x, m = p, n = n)
  I_vect <- rep(1, times = p)
  mu_0 <- as.numeric((t(I_vect) %*% invSS %*% means)/(t(I_vect) %*% 
                                                        invSS %*% I_vect))
  alp_JS_hat <- as.numeric((p + 2)/(p + 2 + n * t(means - 
                                                    mu_0 * I_vect) %*% invSS %*% (means - mu_0 * I_vect)))
  mu_hat_BS <- (1 - alp_JS_hat) * means + alp_JS_hat * mu_0 * 
    I_vect
  list(means = mu_hat_BS, alpha = alp_JS_hat)
}
