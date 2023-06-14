# Purpose: Analysis of weighted random forest results

#  Load packages
library(ggplot2)
library(reshape)

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

# Plot RMSE of all models per number of observations as boxplots
plot <- ggplot(
  results_long,
  aes(x = variable, y = value, fill = variable)
) +
  geom_boxplot() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "RMSE") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_continuous(limits = c(0.5, 2)) +
  scale_x_discrete(labels = c("RF sample cov", "RF nls")) +
  facet_wrap(~n_obs) +
  stat_summary(fun.y = mean, geom = "point", shape = 23, size = 2) +
  ggtitle("Weighted random forest RMSE ratios compared to 1/p weights")
ggsave("results/weighted_rf_rmse_ratios.png", plot)

# Plot RMSE of all models per dataset as boxplots
for (dataset in datasets$datasets) {
  print(
    ggplot(
      results_long[results_long$dataset == dataset, ],
      aes(x = variable, y = value, fill = variable)
    ) +
      geom_boxplot() +
      theme_bw() +
      theme(legend.position = "none") +
      labs(x = "Number of observations", y = "RMSE") +
      facet_wrap(~n_obs, scales = "free_y") +
      stat_summary(fun.y = mean, geom = "point", shape = 23, size = 3) +
      ggtitle(dataset)
  )
}
