# Libraries
library(OpenML)
library(parallel)
source("src/simulations/get_simulation_iid.R")
setOMLConfig(apikey = "c1994bdb7ecb3c6f3c8f3b35f4b47f1f")

# Load names of datasets
datasets <- read.csv("metadata/metadata.csv")

# Set simulation parameters
set.seed(42)
n_cores <- nrow(datasets)
n_trees <- 500 # Ranger default
n_obs <- list(200, 400, 600, 800, 1000, 2000, 3000, 4000, 5000)
n_sim <- 50
# kappas <- list(1, 1.5, 2, 2.5, 100)
kappas <- list(1, 2)

# Run simulations in parallel
mclapply(
  datasets$dataset_id,
  get_simulation_iid,
  mc.cores = n_cores,
  num_trees = n_trees,
  n_obs = n_obs,
  n_sim = n_sim,
  kappas = kappas
)
