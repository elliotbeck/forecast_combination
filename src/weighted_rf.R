# Libraries
library(randomForest)
library(pmlbr)
library(CVXR)
library(reshape)
library(ggplot2)
source("src/rmse.R")
source("src/qis.R")

#  Load benchmark data
data <- fetch_data("225_puma8NH")

# Set simulation parameters
set.seed(1)
n_obs <- list(1000, 2000, 3000, 4000, 5000)
n_sim <- 2

# Set up dataframe to store results
results <- data.frame(
    n_obs = rep(NA, length(n_obs) * n_sim),
    rmse_naive = rep(NA, length(n_obs) * n_sim),
    rmse_lm = rep(NA, length(n_obs) * n_sim),
    rmse_rf = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted = rep(NA, length(n_obs) * n_sim),
    rmse_rf_weighted_shrinkage = rep(NA, length(n_obs) * n_sim)
)

#  Run simulations
k <- 1
for (i in n_obs) {
    for (j in 1:n_sim) {
        # Shuffle rows
        data <- data[sample(1:nrow(data)), ]

        #  Data split
        train_data <- data[1:i, ]
        test_data <- data[(i + 1):(i + 2000), ]

        # Linear regression benchmark
        lm_model <- lm(target ~ ., data = train_data)
        lm_predictions <- predict(lm_model, test_data)

        # Random forest benchmark
        rf_model <- randomForest(target ~ .,
            data = train_data,
            nodesize = 1,
        )

        rf_predictions <- as.data.frame(predict(
            rf_model,
            test_data
        ))

        #  Weighted random forest
        # Calculate weights
        rf_predictions_train <- as.data.frame(predict(
            rf_model,
            train_data,
            predict.all = TRUE
        ))
        rf_residuals <- rf_predictions_train - train_data$target

        # Calculate covariance matrix and mean vector
        mean_vector <- colMeans(rf_residuals)
        cov_matrix <- cov(rf_residuals)
        cov_matrix_shrinkage <- qis(rf_residuals)

        # Calculate weights
        w <- Variable(ncol(rf_residuals))
        objective <- quad_form(w, cov_matrix)
        # objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
        constraints <- list(sum(w) == 1)
        prob <- Problem(Minimize(objective), constraints)
        solution <- solve(prob)
        hist(solution$getValue(w))

        #  Weighted random forest predictions
        rf_predictions_weighted <- as.data.frame(predict(
            rf_model,
            test_data,
            predict.all = TRUE
        ))

        #  Multiply predictions by weights per row
        rf_predictions_weighted <- as.matrix(rf_predictions_weighted) %*%
            as.numeric(solution$getValue(w))

        # Calculate weights
        w <- Variable(ncol(rf_residuals))
        objective <- quad_form(w, cov_matrix_shrinkage)
        # objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
        constraints <- list(sum(w) == 1)
        prob <- Problem(Minimize(objective), constraints)
        solution <- solve(prob)
        hist(solution$getValue(w))

        #  Weighted random forest predictions
        rf_predictions_weighted_shrinkage <- as.data.frame(predict(
            rf_model,
            test_data,
            predict.all = TRUE
        ))

        #  Multiply predictions by weights per row
        rf_predictions_weighted_shrinkage <- as.matrix(rf_predictions_weighted_shrinkage) %*% as.numeric(solution$getValue(w))

        #  Store results
        results[k, ] <- c(
            i,
            rmse(mean(train_data$target), test_data$target), # naive benchmark
            rmse(lm_predictions, test_data$target), # linear regression
            rmse(rf_predictions[, 1], test_data$target), # random forest 1/n weights
            rmse(rf_predictions_weighted, test_data$target), # random forest weighted
            rmse(rf_predictions_weighted_shrinkage, test_data$target) # random forest weighted shrinkage
        )
        k <- k + 1
    }
}

# Get mean results
results_mean <- aggregate(
    results[, 2:6],
    list(results$n_obs),
    mean
)
# Save results
save(results_mean, file = "results/weighted_rf_mean.RData")
save(results, file = "results/weighted_rf_all.RData")

# Plot results df
results_mean_long <- melt(results_mean, id.vars = "Group.1")
ggplot(results_mean_long, aes(x = Group.1, y = value, color = variable)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(
        x = "Number of observations",
        y = "RMSE",
        color = "Model"
    ) +
    scale_x_continuous(breaks = n_obs)
