# Purpose: Analysis of weighted random forest results
#  Load packages
library(ggplot2)
library(reshape)
library(ranger)
library(CVXR)
source("src/utils/qis.R")

# Load names of datasets
datasets <- read.csv("metadata/numerical_regression.csv")
datasets <- datasets[datasets$Remove == 0, ]
datasets <- datasets[datasets$checks_passed == 1, ]
datasets <- datasets[!is.na(datasets$Remove), ]
datasets <- datasets$dataset_name

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
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5

results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:(ncol(results_ratios_rf))
] / results_ratios_rf$mse_rf

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "mse_rf_weighted_1",
    "mse_rf_weighted_shrinkage_1",
    "mse_rf_weighted_1.5",
    "mse_rf_weighted_shrinkage_1.5",
    "mse_rf_weighted_2",
    "mse_rf_weighted_shrinkage_2",
    "mse_rf_weighted_2.5",
    "mse_rf_weighted_shrinkage_2.5",
    "mse_rf_weighted_100",
    "mse_rf_weighted_shrinkage_100"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)

# Calculate various ratios and visualize
result_ratios_rf_long_sample <- results_ratios_rf_long[
  !grepl(results_ratios_rf_long$variable, pattern = "shrinkage"),
]
plot <- ggplot(
  result_ratios_rf_long_sample,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = rep(c("#7CAE00"), 5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0.6, 1.8, .2), limits = c(0.6, 1.8)) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_appendix_sample.pdf", plot)

ratios_rf_long_shrinkage <- results_ratios_rf_long[
  grepl(results_ratios_rf_long$variable, pattern = "shrinkage"),
]
plot <- ggplot(
  ratios_rf_long_shrinkage,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = rep(c("#00BFc4"), 5)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0.6, 1.8, .2), limits = c(0.6, 1.8)) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_appendix_shrinkage.pdf", plot)

results_ratios_rf_long_kappa_2 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_shrinkage_2",
]
plot <- ggplot(
  results_ratios_rf_long_kappa_2,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = "#00BFc4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.7, 1.1), expand = c(0, 0)) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_kappa_2_shrinkage.pdf", plot)

results_ratios_rf_long_kappa_2 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_2",
]
plot <- ggplot(
  results_ratios_rf_long_kappa_2,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = "#7CAE00") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(limits = c(0.7, 1.1), expand = c(0, 0)) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_kappa_2_sample.pdf", plot)

# Calculate and visualize ratios with short-selling vs without short-selling
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf$mse_rf_weighted_2 <- results_ratios_rf$mse_rf_weighted_2 /
  results_ratios_rf$mse_rf_weighted_1
results_ratios_rf$mse_rf_weighted_shrinkage_2 <- results_ratios_rf$mse_rf_weighted_shrinkage_2 /
  results_ratios_rf$mse_rf_weighted_shrinkage_1

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "mse_rf_weighted_2",
    "mse_rf_weighted_shrinkage_2"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)
ratios_kappa_2_sample_1 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_2",
]

plot <- ggplot(
  ratios_kappa_2_sample_1,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00BFc4", "#7CAE00")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.9, 1.05),
    minor_breaks = seq(0.9, 1.05, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_sample_2_sample_1.pdf", plot)

ratios_kappa_2_1_sample <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_shrinkage_2",
]

plot <- ggplot(
  ratios_kappa_2_1_sample,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = "#00BFc4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.9, 1.05),
    minor_breaks = seq(0.9, 1.05, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_rf_long_shrinkage_2_sample_1.pdf", plot)


results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
] / results_ratios_rf$mse_rf_weighted_shrinkage_1

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "mse_rf_weighted_1",
    "mse_rf_weighted_shrinkage_1",
    "mse_rf_weighted_1.5",
    "mse_rf_weighted_shrinkage_1.5",
    "mse_rf_weighted_2",
    "mse_rf_weighted_shrinkage_2",
    "mse_rf_weighted_2.5",
    "mse_rf_weighted_shrinkage_2.5",
    "mse_rf_weighted_100",
    "mse_rf_weighted_shrinkage_100"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)
results_kappa_2_shrinkage_1 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_shrinkage_2",
]

plot <- ggplot(
  results_kappa_2_shrinkage_1,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = "#00BFc4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.9, 1.05),
    minor_breaks = seq(0.9, 1.05, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_kappa_shrinkage_2_shrinkage_1.pdf", plot)

# Calculate and visualize ratios compared to winham et al.
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
] / results_ratios_rf$mse_rf_winham

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "mse_rf_weighted_1",
    "mse_rf_weighted_shrinkage_1",
    "mse_rf_weighted_1.5",
    "mse_rf_weighted_shrinkage_1.5",
    "mse_rf_weighted_2",
    "mse_rf_weighted_shrinkage_2",
    "mse_rf_weighted_2.5",
    "mse_rf_weighted_shrinkage_2.5",
    "mse_rf_weighted_100",
    "mse_rf_weighted_shrinkage_100"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)
results_ratios_rf_winham <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "mse_rf_weighted_shrinkage_2",
]

plot <- ggplot(
  results_ratios_rf_winham,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = "#00BFc4") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.8, 1.1),
    minor_breaks = seq(0.85, 1.05, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_winham.pdf", plot)
