# Purpose: Analysis of weighted random forest results

#  Load packages
library(ggplot2)
library(reshape)
library(pmlbr)
library(ranger)
library(CVXR)
source("src/utils/qis.R")

# Load names of datasets
datasets <- read.table("metadata/datasets.txt", header = TRUE)

#  Load results
load_datasets <- function(dataset) {
  results <- get(load(paste0("results/weighted_rf_unbiased_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results <- lapply(datasets$datasets, load_datasets)
results <- do.call(rbind, results)

# Calculate RMSE ratios
results$rmse_rf_weighted_ratio <- results$rmse_rf_weighted / results$rmse_rf
results$rmse_rf_weighted_shrinkage_ratio <- results$rmse_rf_weighted_shrinkage / results$rmse_rf

#  Convert to long format
results_long <- melt(
  results,
  measure.vars = c(
    "rmse_rf_weighted_ratio",
    "rmse_rf_weighted_shrinkage_ratio"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_long$n_obs <- paste0("n=", results_long$n_obs)
results_long$n_obs <- factor(
  results_long$n_obs,
  levels = paste0("n=", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)

# Plot RMSE of all models per number of observations as boxplots
plot <- ggplot(
  results_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(limits = c(0.5, 2)) +
  scale_x_discrete(labels = c("sample cov", "nls")) +
  facet_wrap(~n_obs) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2)
ggsave("results/weighted_rf_rmse_ratios.pdf", plot)

# Only keep 334_mv dataset
results_long_334_mv <- results_long[results_long$dataset == "344_mv", ]
plot <- ggplot(
  results_long_334_mv,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  # scale_y_continuous(limits = c(0.5, 2)) +
  scale_x_discrete(labels = c("sample cov", "nls")) +
  facet_wrap(~n_obs) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2)
ggsave("results/weighted_rf_rmse_ratios_344_mv.pdf", plot)

#  Plot the weights of the weighted random forest
# Load data
data <- fetch_data(datasets$datasets[1])

# Shuffle rows
data <- data[sample(seq_len(nrow(data))), ]

# Data split
train_data <- data[1:400, ]
test_data <- data[(400 + 1):(400 + 1000), ]

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
weights_sample <- data.frame(weights_sample = solution$getValue(w))

# Calculate weights with shrinkage
w <- Variable(ncol(rf_residuals_train))
objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
constraints <- list(sum(w) == 1)
prob <- Problem(Minimize(objective), constraints)
solution <- solve(prob)
weights_shrinkage <- data.frame(weights_shrinkage = solution$getValue(w))

#  Use ggplot to plot the weights with histogram
weights <- cbind(weights_sample, weights_shrinkage)
weights <- melt(weights)
weights$variable <- factor(weights$variable, levels = c("weights_sample", "weights_shrinkage"))
ggplot(weights, aes(x = value, fill = variable)) +
  geom_histogram(color = "black") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_wrap(~variable, scales = "fixed")
ggsave("results/weighted_rf_weights.pdf")
