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
source("src/simulations/run_iteration.R")

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

  # Limit data to have 10k observations
  data <- data[sample(seq_len(nrow(data)), min(nrow(data), 10000)), ]

  # Run simulations
  results <- mclapply(
    n_obs,
    function(n) {
      mclapply(
        1:n_sim,
        function(i) {
          run_iteration(
            data = data,
            num_trees = num_trees,
            n_obs = n,
            kappas = kappas
          )
        },
        mc.cores = 1
      )
    },
    mc.cores = 1
  )

  # Convert nested list to data frame
  results <- data.frame(do.call(rbind, lapply(results, function(x) do.call(rbind, x))))

  # Get mean results
  results_mean <- aggregate(
    results[, 2:ncol(results)],
    by = list(results$n_obs),
    mean,
  )
  colnames(results_mean)[1] <- "n_obs"

  # Save results
  save(results, file = paste0(
    "results/weighted_rf_", dataset, ".RData"
  ))
  save(results_mean, file = paste0(
    "results/weighted_rf_mean_", dataset, ".RData"
  ))
}
