# Load libraries
library(BVAR)
library(ranger)
library(ggplot2)
source("src/utils/optimizer.R")
source("src/utils/rmse.R")
source("src/utils/qis.R")

# Set target variable
target <- "CPIAUCSL"

# Set kappas
kappas <- list(1, 1.5, 2.5)

# Load data
fred_md_transformed <- fred_transform(fred_md, na.rm = FALSE)

# Transform target variable to log difference
fred_md[, target] <- log(fred_md[, target]) - log(dplyr::lag(fred_md[, target]))
fred_md_transformed[, target] <- fred_md[
    rownames(fred_md) %in% rownames(fred_md_transformed), target
]

# Add lagged target variable
fred_md_transformed <- cbind(
    fred_md_transformed,
    dplyr::lag(fred_md_transformed[, target], 1),
    dplyr::lag(fred_md_transformed[, target], 2),
    dplyr::lag(fred_md_transformed[, target], 3),
    dplyr::lag(fred_md_transformed[, target], 4)
)

# Keep data only after 1980
fred_md_transformed <- fred_md_transformed[
    rownames(fred_md_transformed) >= "1980-01-01",
]

# Keep only non NA columns
fred_md_transformed <- fred_md_transformed[
    ,
    colSums(is.na(fred_md_transformed)) == 0
]

# Print number of remaining predictors
print(paste0(
    "Number of predictors before preprocessing: ",
    ncol(fred_md) - 1
))
print(paste0(
    "Number of predictors after preprocessing: ",
    ncol(fred_md_transformed) - 1
))

# Set predictors and labels
predictors <- fred_md_transformed
labels <- fred_md_transformed[, target]

# Lag predictors
predictors <- predictors[1:(nrow(predictors) - 1), ]
labels <- labels[2:length(labels)]


# Perform rolling window time series cross-validation
window_size <- 200

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
        weight_optimizer,
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
        weight_optimizer,
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
