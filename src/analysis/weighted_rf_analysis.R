# Purpose: Analysis of weighted random forest results

#  Load packages
library(ggplot2)
library(reshape)
library(pmlbr)
library(ranger)
library(CVXR)
source("src/utils/qis.R")

# Load names of datasets
datasets <- summary_stats[summary_stats$task == "regression", ]
datasets <- datasets[
  datasets$n_instances > 6000 & datasets$n_instances < 100000,
]
datasets <- datasets$dataset
datasets <- datasets[datasets != "294_satellite_image"] # Not a regression task

#  Load results
load_datasets <- function(dataset) {
  results <- get(load(paste0("results/weighted_rf_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results <- lapply(datasets, load_datasets)
results <- do.call(rbind, results)

#  Calculate ratios compared to random forest
results_ratios_rf <- results
results_ratios_rf[, 2:(ncol(results_ratios_rf) - 1)] <- results_ratios_rf[
  , 2:(ncol(results_ratios_rf) - 1)
] / results_ratios_rf$rmse_rf

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "rmse_rf_weighted_1",
    "rmse_rf_weighted_shrinkage_1",
    "rmse_rf_weighted_1.5",
    "rmse_rf_weighted_shrinkage_1.5",
    "rmse_rf_weighted_2",
    "rmse_rf_weighted_shrinkage_2",
    "rmse_rf_weighted_2.5",
    "rmse_rf_weighted_shrinkage_2.5",
    "rmse_rf_weighted_100",
    "rmse_rf_weighted_shrinkage_100"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)

# Plot RMSE of all models per number of observations as boxplots
plot <- ggplot(
  results_ratios_rf_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = rep(c("#7CAE00", "#00BFc4"), 5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  ylim(c(0.5, 1.75)) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios.pdf", plot)

# Only keep 334_mv dataset
results_long_334_mv <- results_ratios_rf_long[
  results_ratios_rf_long$dataset == "344_mv",
]
plot <- ggplot(
  results_long_334_mv,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = rep(c("#7CAE00", "#00BFc4"), 5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_344_mv.pdf", plot)

# Compare short/long vs. long only
results_rf_weighted_1 <- results
results_rf_weighted_1[, 2:(ncol(results_rf_weighted_1) - 1)] <-
  results_rf_weighted_1[
    ,
    2:(ncol(results_rf_weighted_1) - 1)
  ] / results_rf_weighted_1$rmse_rf_weighted_1

#  Convert to long format
results_ratios_rf_long <- melt(
  results_rf_weighted_1,
  measure.vars = c(
    "rmse_rf_weighted_shrinkage_1",
    "rmse_rf_weighted_1.5",
    "rmse_rf_weighted_shrinkage_1.5",
    "rmse_rf_weighted_2",
    "rmse_rf_weighted_shrinkage_2",
    "rmse_rf_weighted_2.5",
    "rmse_rf_weighted_shrinkage_2.5"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)

# Plot RMSE of all models per number of observations as boxplots
plot <- ggplot(
  results_ratios_rf_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00BFc4", rep(c("#7CAE00", "#00BFc4"), 3))) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_1.pdf", plot)
