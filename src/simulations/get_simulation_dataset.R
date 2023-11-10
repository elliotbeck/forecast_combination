# Libraries
library(ranger)
library(mlr)
library(tuneRanger)
library(OpenML)
library(reshape)
library(ggplot2)
library(parallel)
source("src/utils/get_mse.R")
source("src/utils/get_qis.R")
source("src/utils/get_winham_benchmark.R")
source("src/utils/get_cesaro_benchmark.R")
source("src/simulations/get_performance.R")

# Run simulations per dataset
get_simulation_dataset_iid <- function(
    dataset,
    num_trees,
    n_obs,
    n_sim,
    kappas) {
  # Load data
  task <- getOMLDataSet(data.name = dataset)
  data <- task$data
  colnames(data)[colnames(data) == task$target.features] <- "target"
  if (dataset == "diamonds") {
    data <- data[, -c(2, 3, 4)]
  }

  # Limit data to have 10k observations
  data <- data[sample(seq_len(nrow(data)), min(nrow(data), 10000)), ]

  # Set up dataframe to store results
  results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    mse_rf = rep(NA, length(n_obs) * n_sim),
    mse_winham = rep(NA, length(n_obs) * n_sim),
    mse_cesaro = rep(NA, length(n_obs) * n_sim)
  )

  results <- cbind(
    results,
    data.frame(
      matrix(
        rep(NA, length(n_obs) * n_sim * (length(kappas) * 2)),
        nrow = length(n_obs) * n_sim,
        ncol = length(kappas) * 2,
        dimnames = list(
          NULL,
          c(
            paste0("mse_rf_weighted_", kappas),
            paste0("mse_rf_weighted_shrinkage_", kappas)
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
      test_data <- data[(i + 1):nrow(data), ]

      # Standardizing target variable, not allowing for leakage
      norm_param <- list(
        mean = mean(train_data$target),
        sd = sd(train_data$target)
      )
      train_data$target <- (train_data$target - norm_param$mean) / norm_param$sd

      rf_model <- ranger(
        target ~ .,
        data = train_data,
        num.trees = num_trees,
        mtry = round((ncol(train_data) - 1) / 3),
        replace = TRUE,
        keep.inbag = TRUE
      )

      rf_predictions <- predict(rf_model, test_data)$predictions
      rf_predictions <- (rf_predictions * norm_param$sd) + norm_param$mean

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

      # Get mse for all kappas and sample cov matrix
      results_sample <- sapply(
        kappas,
        get_performance,
        mean_vector = mean_vector,
        cov_matrix = cov_matrix,
        predictions_test = rf_predictions_test_all,
        labels_test = test_data$target,
        mean = norm_param$mean,
        sd = norm_param$sd
      )

      # Get mse for all kappas and nls cov matrix
      results_nls <- sapply(
        kappas,
        get_performance,
        mean_vector = mean_vector,
        cov_matrix = cov_matrix_shrinkage,
        predictions_test = rf_predictions_test_all,
        labels_test = test_data$target,
        mean = norm_param$mean,
        sd = norm_param$sd
      )

      # Bind all results together
      mse_bkw <- c(results_sample, results_nls)

      # Get winham et al. benchmark
      mse_winham <- winham(
        train_data = train_data,
        test_data = test_data,
        rf_model = rf_model,
        norm_param = norm_param,
        rf_predictions_train_all = rf_predictions_train_all,
        rf_predictions_test_all = rf_predictions_test_all
      )

      # Get cesaro benchmark
      mse_cesaro <- cesaro(
        train_data = train_data,
        test_data = test_data,
        rf_model = rf_model,
        norm_param = norm_param,
        rf_predictions_train_all = rf_predictions_train_all,
        rf_predictions_test_all = rf_predictions_test_all
      )

      #  Store results
      results[k, ] <- c(
        i,
        mse(rf_predictions, test_data$target),
        mse_winham,
        mse_cesaro,
        mse_bkw
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
