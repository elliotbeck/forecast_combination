# Libraries
library(ranger)
library(pmlbr)
library(CVXR)
library(reshape)
library(ggplot2)
source("src/utils/rmse.R")
source("src/utils/qis.R")

# Load names of datasets
datasets <- read.table("metadata/datasets.txt", header = TRUE)

# Set simulation parameters
set.seed(1)
n_obs <- list(200, 400, 600, 1000, 2000, 3000, 4000, 5000)
n_sim <- 50

# Run simulations per dataset
for (dataset in datasets$datasets) {
  # Load data
  data <- fetch_data(dataset)

  # Set up dataframe to store results
  results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    rmse_naive = rep(NA, length(n_obs) * n_sim),
    rmse_lm = rep(NA, length(n_obs) * n_sim),
    rmse_rf = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage = rep(NA, length(n_obs) * n_sim)
  )

  # Run simulations
  k <- 1
  for (i in n_obs) {
    for (j in 1:n_sim) {
      # Shuffle rows
      data <- data[sample(seq_len(nrow(data))), ]

      # Data split
      train_data <- data[1:i, ]
      test_data <- data[(i + 1):(i + 1000), ]

      # Linear regression benchmark
      lm_model <- lm(target ~ ., data = train_data)
      lm_predictions <- predict(lm_model, test_data)

      # Random forest benchmark with ranger
      rf_model <- ranger(target ~ ., data = train_data)
      rf_predictions <- predict(rf_model, test_data)$predictions

      #  Weighted random forest
      rf_predictions_train_all <- predict(
        rf_model, train_data,
        predict.all = TRUE
      )$predictions
      rf_residuals_train <- train_data$target - rf_predictions_train_all

      # Calculate covariance matrix and mean vector of residuals
      mean_vector <- colMeans(rf_residuals_train)
      cov_matrix <- cov(rf_residuals_train)
      cov_matrix_shrinkage <- qis(rf_residuals_train)

      # Calculate weights without shrinkage
      w <- Variable(ncol(rf_residuals_train))
      objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix)
      constraints <- list(sum(w) == 1)
      prob <- Problem(Minimize(objective), constraints)
      solution <- solve(prob)

      #  Weighted random forest predictions
      rf_predictions_test_all <- predict(rf_model, test_data, predict.all = TRUE)$predictions

      # #  Multiply predictions by weights per row
      rf_predictions_weighted <- as.matrix(rf_predictions_test_all) %*% solution$getValue(w)

      # Calculate weights with shrinkage
      w <- Variable(ncol(rf_residuals_train))
      objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
      constraints <- list(sum(w) == 1)
      prob <- Problem(Minimize(objective), constraints)
      solution <- solve(prob)

      # Shrinkage weighted random forest predictions
      rf_predictions_weighted_nls <- as.matrix(rf_predictions_test_all) %*%
        solution$getValue(w)

      #  Store results
      results[k, ] <- c(
        i,
        rmse(mean(train_data$target), test_data$target),
        rmse(lm_predictions, test_data$target),
        rmse(rf_predictions, test_data$target),
        rmse(rf_predictions_weighted, test_data$target),
        rmse(rf_predictions_weighted_nls, test_data$target)
      )
      k <- k + 1
      print(k)
    }
  }

  # Get mean results
  results_mean <- aggregate(
    results[, 2:6],
    list(results$n_obs),
    mean
  )
  # Save results
  save(results, file = paste0("results/weighted_rf_unbiased_", dataset, ".RData"))
  save(results_mean, file = paste0("results/weighted_rf_mean_unbiased_", dataset, ".RData"))

  # Convert results to long format
  results_mean_long <- melt(results_mean, id.vars = "Group.1")
  results_long <- melt(results, id.vars = "n_obs")

  # Plot mean results df
  plot <- ggplot(results_mean_long, aes(x = Group.1, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(
      x = "Number of observations",
      y = "RMSE",
      color = "Model"
    )
  ggsave(paste0("results/weighted_rf_unbiased_", dataset, ".png"), plot)
}
