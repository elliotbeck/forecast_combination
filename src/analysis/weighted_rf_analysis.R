# Purpose: Analysis of weighted random forest results
#  Load packages
library(ggplot2)
library(reshape)
library(xtable)
source("src/utils/print_bold.R")

# Load names of datasets
datasets <- read.csv("metadata/numerical_regression.csv")

#  Load results
load_datasets <- function(dataset) {
  results <- get(load(paste0("results/weighted_rf_mean_", dataset, ".RData")))
  results$dataset <- dataset
  return(results)
}
results <- lapply(datasets$dataset_id, load_datasets)
results <- do.call(rbind, results)

#  Calculate ratios compared to random forest
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5

results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:(ncol(results_ratios_rf))
] / results_ratios_rf$rf

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "sample_1",
    "nls_1",
    "sample_1.5",
    "nls_1.5",
    "sample_2",
    "nls_2",
    "sample_2.5",
    "nls_2.5",
    "sample_100",
    "nls_100"
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
  !grepl(results_ratios_rf_long$variable, pattern = "nls"),
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
ggsave("results/weighted_rf_rmse_ratios_appendix_sample.eps", plot)

ratios_rf_long_shrinkage <- results_ratios_rf_long[
  grepl(results_ratios_rf_long$variable, pattern = "nls"),
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
ggsave("results/weighted_rf_rmse_ratios_appendix_shrinkage.eps", plot)

results_ratios_rf_long_kappa_2 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "nls_2",
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
ggsave("results/weighted_rf_rmse_ratios_kappa_2_shrinkage.eps", plot)

results_ratios_rf_long_kappa_2 <- results_ratios_rf_long[
  results_ratios_rf_long$variable == "sample_2",
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
ggsave("results/weighted_rf_rmse_ratios_kappa_2_sample.eps", plot)

# Calculate and visualize ratios with short-selling vs without short-selling
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf$sample_2 <- results_ratios_rf$sample_2 /
  results_ratios_rf$sample_1
results_ratios_rf$nls_2 <- results_ratios_rf$nls_2 /
  results_ratios_rf$nls_1

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "sample_2",
    "nls_2"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)
plot <- ggplot(
  results_ratios_rf_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#7CAE00", "#00BFc4")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.9, 1.05),
    minor_breaks = seq(0.9, 1.05, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_short_selling.eps", plot)

# Calculate and visualize ratios compared to sample 1
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf$mse_rf_weighted_2 <- results_ratios_rf$sample_2 /
  results_ratios_rf$sample_1
results_ratios_rf$mse_rf_weighted_shrinkage_2 <- results_ratios_rf$nls_2 /
  results_ratios_rf$sample_1

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
ggsave("results/weighted_rf_rmse_ratios_sample_2_sample_1.eps", plot)

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
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_shrinkage_2_sample_1.eps", plot)

# Calculate and visualize ratios compared to winham et al. and pham et al.
results_ratios_rf <- results
results_ratios_rf <- aggregate(. ~ n_obs + dataset, results_ratios_rf, mean)
results_ratios_rf[, 3:ncol(results_ratios_rf)] <- results_ratios_rf[
  , 3:ncol(results_ratios_rf)
]^0.5
results_ratios_rf$winham <- results_ratios_rf$nls_2 /
  results_ratios_rf$wrf
results_ratios_rf$cesaro <- results_ratios_rf$nls_2 /
  results_ratios_rf$crf

#  Convert to long format
results_ratios_rf_long <- melt(
  results_ratios_rf,
  measure.vars = c(
    "winham",
    "cesaro"
  ),
  id.vars = c("dataset", "n_obs"),
)

# Change names of number of observations
results_ratios_rf_long$n_obs <- paste0("n = ", results_ratios_rf_long$n_obs)
results_ratios_rf_long$n_obs <- factor(
  results_ratios_rf_long$n_obs,
  levels = paste0("n = ", c(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000))
)
plot <- ggplot(
  results_ratios_rf_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_blank()) +
  scale_fill_manual(values = c("#00B9E3", "#619CFF")) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  scale_y_continuous(
    limits = c(0.8, 1.1),
    minor_breaks = seq(0.85, 1.1, 0.025)
  ) +
  facet_wrap(~n_obs, nrow = 1, strip.position = "bottom") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "red")
ggsave("results/weighted_rf_rmse_ratios_benchmarks.eps", plot)

#  Get RMSE tables for each number of observations
results_rmse <- results
results_rmse <- aggregate(. ~ n_obs + dataset, results_rmse, mean)
results_rmse[, 3:ncol(results_rmse)] <- results_rmse[
  , 3:ncol(results_rmse)
]^0.5

results_rmse$dataset_name <- datasets$dataset_name[match(
  results_rmse$dataset,
  datasets$dataset_id
)]

for (n_obs in unique(results_rmse$n_obs)) {
  results_rmse_subset <- results_rmse[results_rmse$n_obs == n_obs, ]
  # order alphabetically
  results_rmse_subset <- results_rmse_subset[
    order(results_rmse_subset$dataset_name),
  ]
  results_rmse_subset <- results_rmse_subset[
    c("dataset_name", "rf", "nls_2", "wrf", "crf")
  ]
  colnames(results_rmse_subset) <- c("Name", "RF", "HRF", "WRF", "CRF")
  results_rmse_subset[
    results_rmse_subset$Name %in% c("Ailerons", "elevators"), 2:ncol(results_rmse_subset)
  ] <- results_rmse_subset[
    results_rmse_subset$Name %in% c("Ailerons", "elevators"), 2:ncol(results_rmse_subset)
  ] * 1000
  table <- xtable(
    digits = 5,
    results_rmse_subset,
    caption = paste0(
      "Results for all data sets with n = ",
      n_obs,
      ". Values for the Ailerons and elevators data sets are scaled by a factor of 1,000 for better readability."
    )
  )
  table <- printbold(
    table,
    each = "row",
    max = FALSE,
    file = paste0("results/tables/weighted_rf_rmse_", n_obs, ".tex"),
    include.rownames = FALSE
  )
  print(table, type = "latex", file = paste0("results/tables/weighted_rf_rmse_", n_obs, ".tex"))
}
