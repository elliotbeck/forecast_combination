# Libraries
library(ranger)
library(mlr)
library(tuneRanger)
library(OpenML)
library(reshape)
library(ggplot2)
library(parallel)
source("src/utils/get_mse.R")
source("src/cov_estimators/get_cov_qis.R")
source("src/utils/get_winham_benchmark.R")
source("src/utils/get_cesaro_benchmark.R")
source("src/utils/get_performance.R")

# Run simulations per dataset
get_simulation_iid <- function(
    dataset,
    num_trees,
    n_obs,
    n_sim,
    kappas) {
  # Load data
  task <- getOMLDataSet(data.id = dataset)
  data <- task$data
  colnames(data)[colnames(data) == task$target.features] <- "target"
  if (dataset == "diamonds") {
    data <- data[, -c(2, 3, 4)]
  }

  # Load hyperparameters
  hyperparams <- read.csv("metadata/metadata.csv")
  hyperparams <- hyperparams[hyperparams$dataset_id == dataset, ]

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
    print(paste0("Dataset: ", dataset, ", n_obs: ", i))
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

      # Train default random forest
      rf_model <- ranger(
        target ~ .,
        data = train_data,
        num.trees = num_trees,
        mtry = floor((ncol(train_data) - 1) / 3),
        replace = TRUE,
        keep.inbag = TRUE
      )

      # # Fit random forest with tuned mtry as benchmark
      # rf_model_mtry <- tuneRanger::tuneMtryFast(
      #   target ~ .,
      #   data = train_data,
      #   num.treesTry = num_trees,
      #   mtryStart = round((ncol(train_data) - 1) / 3),
      #   doBest = TRUE,
      #   plot = FALSE,
      #   trace = FALSE
      # )

      # # Get tuned mtry
      # task <- mlr::makeRegrTask(
      #   data = train_data,
      #   target = "target"
      # )

      # rf_model_mtry <- suppressWarnings((tuneRanger::tuneRanger(
      #   task = task,
      #   num.trees = num_trees,
      #   show.info = FALSE
      # )))

      mtry <- if (hyperparams$model__max_features == "sqrt") {
        floor(sqrt(ncol(train_data) - 1))
      } else if (hyperparams$model__max_features == "log2") {
        floor(log2(ncol(train_data) - 1))
      } else {
        floor((ncol(train_data) - 1) * as.numeric(hyperparams$model__max_features))
      }

      rf_model_tuned <- ranger(
        target ~ .,
        data = train_data,
        num.trees = num_trees,
        mtry = mtry,
        min.node.size = hyperparams$model__min_samples_split,
        replace = if (hyperparams$model__bootstrap == 1) TRUE else FALSE
      )

      # Tuned random forest predictions on test data
      rf_predictions <- predict(rf_model_tuned, test_data)$predictions
      rf_predictions <- (rf_predictions * norm_param$sd) + norm_param$mean

      # Random forest predictions and residuals on train data
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

      # # Mtry hrf
      # mtry_sequence <- runif(num_trees, 0.2, 0.8)
      # min_node_size_sequence <- sample(c(1, 3, 5, 7, 9, 11), num_trees, replace = TRUE)

      # ranger_mtry <- function(mtry, min_node_size, data) {
      #   fit <- ranger(target ~ .,
      #     data = data,
      #     mtry = mtry,
      #     min.node.size = min_node_size,
      #     num.trees = 1
      #   )
      #   return(fit)
      # }
      # predict_mtry <- function(fit, data) {
      #   rf_predictions_all <- predict(
      #     fit,
      #     data,
      #     predict.all = TRUE
      #   )$predictions
      #   return(rf_predictions_all)
      # }
      # fit_ranger_mtry <- mcmapply(
      #   FUN = ranger_mtry,
      #   mtry = mtry_sequence,
      #   min_node_size = min_node_size_sequence,
      #   MoreArgs = list(data = train_data),
      #   SIMPLIFY = FALSE,
      #   mc.cores = 32
      # )
      # rf_predictions_train_all <- mclapply(
      #   fit_ranger_mtry,
      #   predict_mtry,
      #   data = train_data,
      #   mc.cores = 32
      # )
      # rf_predictions_train_all <- do.call(cbind, rf_predictions_train_all)
      # rf_residuals_train <- train_data[, "target"] - rf_predictions_train_all

      # rf_predictions_test_all <- mclapply(
      #   fit_ranger_mtry,
      #   predict_mtry,
      #   data = test_data,
      #   mc.cores = 32
      # )
      # rf_predictions_test_all <- do.call(cbind, rf_predictions_test_all)

      # Calculate covariance matrix and mean vector of residuals
      mean_vector <- colMeans(rf_residuals_train)
      cov_matrix <- cov(rf_residuals_train)
      cov_matrix_shrinkage <- get_cov_qis(rf_residuals_train)

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
