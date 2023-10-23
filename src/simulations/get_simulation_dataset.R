# Libraries
library(ranger)
library(OpenML)
library(reshape)
library(ggplot2)
library(parallel)
library(xgboost)
source("src/utils/mse.R")
source("src/utils/qis.R")
source("src/utils/winham.R")
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

  # Set up dataframe to store results
  results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    mse_rf = rep(NA, length(n_obs) * n_sim),
    mse_xgb = rep(NA, length(n_obs) * n_sim)
  )

  results <- cbind(
    results,
    data.frame(
      matrix(
        rep(NA, length(n_obs) * n_sim * (length(kappas) * 2 + 1)),
        nrow = length(n_obs) * n_sim,
        ncol = length(kappas) * 2 + 1,
        dimnames = list(
          NULL,
          c(
            paste0("mse_rf_weighted_", kappas),
            paste0("mse_rf_weighted_shrinkage_", kappas),
            "mse_rf_winham"
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

      # Random forest benchmark with ranger
      rf_model <- ranger(
        target ~ .,
        data = train_data,
        num.trees = num_trees,
        keep.inbag = TRUE
      )
      rf_predictions <- predict(rf_model, test_data)$predictions
      rf_predictions <- (rf_predictions * norm_param$sd) + norm_param$mean

      # XGBoost benchmark with xgboost
      xgb_model <- xgboost(
        data = as.matrix(subset(train_data, select = -target)),
        label = train_data$target,
        nrounds = num_trees,
        nthread = 1,
        verbose = 0
      )
      xgb_predictions <- predict(
        xgb_model,
        as.matrix(subset(test_data, select = -target))
      )
      xgb_predictions <- (xgb_predictions * norm_param$sd) + norm_param$mean

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

      #  Store results
      results[k, ] <- c(
        i,
        mse(rf_predictions, test_data$target),
        mse(xgb_predictions, test_data$target),
        mse_bkw,
        mse_winham
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
