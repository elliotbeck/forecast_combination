# Libraries
library(pmlbr)
library(parallel)
source("src/simulations/get_simulation_dataset.R")

# Load names of datasets
datasets <- summary_stats[summary_stats$task == "regression", ]
datasets <- datasets[
  datasets$n_instances >= 6000 & datasets$n_instances <= 100000,
]
datasets <- datasets$dataset
datasets <- datasets[datasets != "294_satellite_image"] # Not a regression task
datasets <- datasets[datasets != "503_wind"] # Time series, not iid

# Set simulation parameters
set.seed(42)
num_trees <- 500 # Ranger default
n_obs <- list(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000)
n_sim <- 100
kappas <- list(1, 1.5, 2, 2.5, 100)

# Run simulations in parallel
mclapply(
  datasets,
  get_simulation_dataset_iid,
  mc.cores = length(datasets),
  num_trees = num_trees,
  n_obs = n_obs,
  n_sim = n_sim,
  kappas = kappas
)
