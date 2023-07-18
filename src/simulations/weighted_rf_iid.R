# Libraries
library(ranger)
library(pmlbr)
library(reshape)
library(ggplot2)
source("src/utils/rmse.R")
source("src/utils/qis.R")
source("src/utils/get_performance.R")

# Load names of datasets
datasets <- summary_stats[summary_stats$task == "regression", ]
datasets <- datasets[
  datasets$n_instances > 6000 & datasets$n_instances < 100000,
]
datasets <- datasets$dataset
datasets <- datasets[datasets != "294_satellite_image"] # Not a regression task

# Set simulation parameters
set.seed(1)
num_trees <- 500 # Ranger default
n_obs <- list(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000)
n_sim <- 50
kappas <- list(1, 1.5, 2, 2.5, 100)

# Run simulations per dataset
for (dataset in datasets) {
  # Load data
  data <- fetch_data(dataset)

  # Set up dataframe to store results
  results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    rmse_rf = rep(NA, length(n_obs) * n_sim)
  )
  results <- cbind(
    results,
    data.frame(
      matrix(
        rep(NA, length(n_obs) * n_sim * length(kappas)),
        nrow = length(n_obs) * n_sim,
        ncol = length(kappas) * 2,
        dimnames = list(
          NULL,
          c(
            paste0("rmse_rf_weighted_", kappas),
            paste0("rmse_rf_weighted_shrinkage_", kappas)
          )
        )
      )
    )
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

      # Normalize target variable, not allowing leakage
      norm_param <- list(
        mean = mean(train_data$target),
        sd = sd(train_data$target)
      )
      train_data$target <- (train_data$target - norm_param$mean) / norm_param$sd
      test_data$target <- (test_data$target - norm_param$mean) / norm_param$sd

      # Random forest benchmark with ranger
      rf_model <- ranger(target ~ ., data = train_data, num.trees = num_trees)
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
        get_performance,
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
        get_performance,
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
        results_all_test
      )
      k <- k + 1
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
}
