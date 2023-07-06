# Libraries
library(ranger)
library(pmlbr)
library(reshape)
library(ggplot2)
source("src/utils/rmse.R")
source("src/utils/qis.R")
source("src/utils/optimizer.R")

# Load names of datasets
datasets <- read.table("metadata/datasets.txt", header = TRUE)

# Set simulation parameters
set.seed(1)
n_obs <- list(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000)
n_sim <- 50
kappas <- list(1, 1.5, 2, 2.5)

# Run simulations per dataset
for (dataset in datasets$datasets) {
  # Load data
  data <- fetch_data(dataset)

  # Set up dataframe to store results
  results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    rmse_rf = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_1 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_1_5 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_2 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_2.5 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage_1 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage_1_5 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage_2 = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage_2.5 = rep(NA, length(n_obs) * n_sim)
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

      # Random forest benchmark with ranger
      rf_model <- ranger(target ~ ., data = train_data)
      rf_predictions <- predict(rf_model, test_data)$predictions

      # Reandom forest predictions and residuals on train data
      rf_predictions_train_all <- predict(
        rf_model,
        train_data,
        predict.all = TRUE
      )$predictions
      rf_residuals_train <- train_data$target - rf_predictions_train_all

      # Random forest predictions on test data
      rf_predictions_test_all <- predict(
        rf_model,
        test_data,
        predict.all = TRUE
      )$predictions

      # Calculate covariance matrix and mean vector of residuals
      mean_vector <- colMeans(rf_residuals_train)
      cov_matrix <- cov(rf_residuals_train)
      cov_matrix_shrinkage <- qis(rf_residuals_train)

      # Get rmse for all kappas and sample cov matrix
      results_sample <- sapply(
        kappas,
        weight_optimizer,
        mean_vector = mean_vector,
        cov_matrix = cov_matrix,
        predictions_train = rf_predictions_train_all,
        predictions_test = rf_predictions_test_all,
        labels_train = train_data$target,
        labels_test = test_data$target
      )

      # Get rmse for all kappas and nls cov matrix
      results_nls <- sapply(
        kappas,
        weight_optimizer,
        mean_vector = mean_vector,
        cov_matrix = cov_matrix_shrinkage,
        predictions_train = rf_predictions_train_all,
        predictions_test = rf_predictions_test_all,
        labels_train = train_data$target,
        labels_test = test_data$target
      )

      # Bind all results together
      results_all_train <- c(results_sample[1, ], results_nls[1, ])
      results_all_test <- c(results_sample[2, ], results_nls[2, ])

      #  Store results
      results[k, ] <- c(
        i,
        rmse(rf_predictions, test_data$target),
        results_all_test[1],
        results_all_test[2],
        results_all_test[3],
        results_all_test[4],
        results_all_test[5],
        results_all_test[6],
        results_all_test[7],
        results_all_test[8]
      )
      k <- k + 1
      print(k)
    }
  }

  # Get mean results
  results_mean <- aggregate(
    results[, 2:ncol(results)],
    list(results$n_obs),
    mean
  )
  # Save results
  save(results, file = paste0(
    "results/weighted_rf_", dataset, ".RData"
  ))
  save(results_mean, file = paste0(
    "results/weighted_rf_mean_", dataset, ".RData"
  ))

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
  ggsave(paste0("results/weighted_rf_", dataset, ".png"), plot)
}
