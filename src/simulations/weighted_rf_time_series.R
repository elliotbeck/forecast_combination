# Load libraries
library(pmlbr)
library(ranger)
library(ggplot2)
source("src/simulations/get_performance.R")
source("src/utils/mse.R")
source("src/utils/qis.R")

# Set target variable
dataset <- "503_wind"

# Load data
data <- fetch_data(dataset)

# Remove date and time variables
data <- data[, -c(1, 2, 3)]

# Set predictors and labels
predictors <- data
predictors[, "target"] <- NULL
labels <- data[, "target"]

# Perform rolling window time series cross-validation
window_size <- 200
kappas <- list(1, 1.5, 2.5)

# Prepare results data frame
results <- data.frame(
  mse_rf = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_sample_1 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_sample_1_5 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_sample_2_5 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_nls_1 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_nls_1_5 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_nls_2_5 = rep(NA, length(window_size:(nrow(predictors) - 1))),
  mse_rw = rep(NA, length(window_size:(nrow(predictors) - 1)))
)
k <- 1
for (i in window_size:(nrow(predictors) - 1)) {
  # Get train and test data
  predictors_train <- predictors[(i - window_size + 1):i, ]
  predictors_test <- predictors[i + 1, ]
  labels_train <- labels[(i - window_size + 1):i]
  labels_test <- labels[i + 1]

  # Normalize target variable, not allowing leakage
  norm_param <- list(
    mean = mean(labels_train),
    sd = sd(labels_train)
  )
  labels_train <- (labels_train - norm_param$mean) / norm_param$sd
  labels_test <- (labels_test - norm_param$mean) / norm_param$sd

  # Train random forest
  rf <- ranger(x = predictors_train, y = labels_train)
  rf_predictions <- predict(rf, predictors_test)$predictions

  # Random forest predictions and residuals on train data
  rf_predictions_train_all <- predict(
    rf,
    predictors_train,
    predict.all = TRUE
  )$predictions
  rf_residuals_train <- labels_train - rf_predictions_train_all

  # Random forest predictions on test data
  rf_predictions_test_all <- predict(
    rf,
    predictors_test,
    predict.all = TRUE
  )$predictions

  # Get mean vector and covariance matrix
  mean_vector <- colMeans(rf_residuals_train)
  cov_matrix <- cov(rf_residuals_train)
  cov_matrix_shrinkage <- qis(rf_residuals_train)

  # Get rmse for all kappas and nls cov matrix
  results_sample <- sapply(
    kappas,
    get_performance,
    mean_vector = mean_vector,
    cov_matrix = cov_matrix,
    predictions_train = rf_predictions_train_all,
    predictions_test = rf_predictions_test_all,
    labels_train = labels_train,
    labels_test = labels_test
  )

  # Get rmse for all kappas and nls cov matrix
  results_nls <- sapply(
    kappas,
    get_performance,
    mean_vector = mean_vector,
    cov_matrix = cov_matrix_shrinkage,
    predictions_train = rf_predictions_train_all,
    predictions_test = rf_predictions_test_all,
    labels_train = labels_train,
    labels_test = labels_test
  )

  # Bind all results together
  results_all <- c(results_sample[2, ], results_nls[2, ])

  # Store results
  results[k, ] <- c(
    rmse(rf_predictions, labels_test),
    results_all[1],
    results_all[2],
    results_all[3],
    results_all[4],
    results_all[5],
    results_all[6],
    rmse(labels_train[length(labels_train)], labels_test)
  )
  print(k)
  k <- k + 1
}

# Get mean of results
results_mean <- colMeans(results, na.rm = TRUE)
print(results_mean)

# Plot cummulative mse
results_cum_sum <- cumsum(results)
results_cum_sum$date <- as.Date(
  rownames(fred_md_transformed)[(window_size + 2):(nrow(fred_md_transformed))]
)

# Transform to long
results_cum_sum_long <- reshape2::melt(
  results_cum_sum,
  id.vars = "date"
)

# Plot
ggplot(results_cum_sum_long, aes(x = date, y = value, color = variable)) +
  geom_line() +
  theme_minimal() +
  labs(
    x = "Observation",
    y = "Cumulative MSE",
    color = "Method"
  ) +
  theme(
    legend.position = "bottom"
  )
