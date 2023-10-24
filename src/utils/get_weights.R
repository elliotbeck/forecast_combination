# Libraries
library(CVXR)

get_weights <- function(mean_vector, cov_matrix, kappa) {
  w <- Variable(ncol(cov_matrix))
  objective <- square((t(w) %*% mean_vector)) + quad_form(w, cov_matrix)
  constraints <- list(
    sum(w) == 1,
    sum(abs(w)) <= kappa
  )
  prob <- Problem(Minimize(objective), constraints)
  solution <- solve(prob, num_iter = 100000)
  return(solution$getValue(w))
}
