# Library
library(dplyr)

# Load names of datasets
datasets <- read.csv("metadata/numerical_regression.csv")
datasets <- datasets[datasets$Remove == 0, ]
datasets <- datasets[datasets$checks_passed == 1, ]
datasets <- datasets[!is.na(datasets$Remove), ]
datasets <- datasets$dataset_name

# Read in benchmark results
hyperparams <- read.csv("metadata/benchmark_total.csv")
metadata <- read.csv("metadata/regression_table.csv")

# Filter for RandomForest and datasets
hyperparams_subset <- hyperparams[
  (hyperparams$data__keyword %in% datasets) &
    (hyperparams$model_name == "RandomForest") &
    (hyperparams$model__criterion == "squared_error") &
    (hyperparams$regression == 1),
]

#  Only keep relevant columns
hyperparams_subset <- hyperparams_subset[
  ,
  c(
    "data__keyword",
    "mean_test_score",
    "model__criterion",
    "model__min_samples_split",
    "model__max_features",
    "model__bootstrap"
  )
]

# Only keep complete cases
hyperparams_subset <- hyperparams_subset[complete.cases(hyperparams_subset), ]

# Order performance by dataset and pick best (only first) hyperparameters
hyperparams_best <- hyperparams_subset %>%
  group_by(data__keyword) %>%
  filter(mean_test_score == max(mean_test_score)) %>%
  slice(1)

# Append dataset_id to hyperparameters (only this column)
hyperparams_best <- merge(
  hyperparams_best,
  metadata[, c("dataset_name", "dataset_id")],
  by.x = "data__keyword",
  by.y = "dataset_name",
  all.x = TRUE
)

# Save hyperparameters
write.csv(hyperparams_best, "metadata/metadata.csv", row.names = FALSE)
