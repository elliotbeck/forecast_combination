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
  results <- get(load(paste0("results/weighted_rf_unbiased_long_only_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results <- lapply(datasets$datasets, load_datasets)
results <- do.call(rbind, results)

# Calculate RMSE ratios
results$rmse_rf_weighted_ratio_1 <- results$rmse_rf_weighted_1 / results$rmse_rf
results$rmse_rf_weighted_ratio_1_6 <- results$rmse_rf_weighted_1_6 / results$rmse_rf
results$rmse_rf_weighted_ratio_100 <- results$rmse_rf_weighted_100 / results$rmse_rf
results$rmse_rf_weighted_shrinkage_ratio_1 <- results$rmse_rf_weighted_shrinkage_1 / results$rmse_rf
results$rmse_rf_weighted_shrinkage_ratio_1_6 <- results$rmse_rf_weighted_shrinkage_1_6 /
  results$rmse_rf
results$rmse_rf_weighted_shrinkage_ratio_100 <- results$rmse_rf_weighted_shrinkage_100 /
  results$rmse_rf

#  Convert to long format
results_long <- melt(
  results,
  measure.vars = c(
    "rmse_rf_weighted_ratio_1",
    "rmse_rf_weighted_shrinkage_ratio_1",
    "rmse_rf_weighted_ratio_1_6",
    "rmse_rf_weighted_shrinkage_ratio_1_6",
    "rmse_rf_weighted_ratio_100",
    "rmse_rf_weighted_shrinkage_ratio_100"
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
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = rep(c("#7CAE00", "#00BFc4", 3))) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.4, 2)) +
  scale_x_discrete(labels = c(rep("sample cov", 3), rep("nls", 3))) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "red")
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
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#7CAE00", "#00BFc4")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  # scale_y_continuous(limits = c(0.4, 2)) +
  scale_x_discrete(labels = c("sample cov", "nls")) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "red")
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
  scale_fill_manual(values = c("#7CAE00", "#00BFc4")) +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_wrap(~variable, scales = "fixed")
ggsave("results/weighted_rf_weights_long_only_.pdf")

# Compare short/long vs. only long
# Load results long only
load_datasets_long <- function(dataset) {
  results <- get(load(paste0("results/weighted_rf_unbiased_long_only_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results_long <- lapply(datasets$datasets, load_datasets_long)
results_long <- do.call(rbind, results_long)
# Load results short
load_datasets_short <- function(dataset) {
  results <- get(load(paste0("results/weighted_rf_unbiased_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results_short <- lapply(datasets$datasets, load_datasets_short)
results_short <- do.call(rbind, results_short)

# Calculate mean RMSEs
results_long_mean <- aggregate(
  results_long[, c(
    "rmse_naive", "rmse_lm", "rmse_rf",
    "rmse_rf_weighted", "rmse_rf_weighted_shrinkage"
  )],
  by = list(results_long$n_obs),
  FUN = mean
)
results_long_mean$method <- "long"
results_short_mean <- aggregate(
  results_short[, c(
    "rmse_naive", "rmse_lm", "rmse_rf",
    "rmse_rf_weighted", "rmse_rf_weighted_shrinkage"
  )],
  by = list(results_short$n_obs),
  FUN = mean
)
results_short_mean$method <- "short"
