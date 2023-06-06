# Libraries
library(randomForest)
library(pmlbr)
library(CVXR)
source("code/rmse.R")
source("code/qis.R")

#  Load benchmark data
data <- fetch_data("225_puma8NH")
dim(data)

# Shuffle rows
data <- data[sample(1:nrow(data)), ]

#  Data split
train_data <- data[1:2000, ]
test_data <- data[6001:7500, ]

# Linear regression benchmark
lm_model <- lm(target ~ ., data = train_data)
lm_predictions <- predict(lm_model, test_data)

# Random forest benchmark
rf_model <- randomForest(target ~ .,
    data = train_data,
    nodesize = 1,
)

rf_predictions <- as.data.frame(predict(
    rf_model,
    test_data
    # predict.all = TRUE
))

#  Weighted random forest
# Calculate weights
rf_predictions_train <- as.data.frame(predict(
    rf_model,
    train_data,
    predict.all = TRUE
))
rf_residuals <- rf_predictions_train - train_data$target

# Calculate covariance matrix and mean vector
mean_vector <- colMeans(rf_residuals)
cov_matrix <- cov(rf_residuals)
cov_matrix_shrinkage <- qis(rf_residuals)

# Calculate weights
w <- Variable(ncol(rf_residuals))
objective <- quad_form(w, cov_matrix)
# objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
constraints <- list(sum(w) == 1)
prob <- Problem(Minimize(objective), constraints)
solution <- solve(prob)
hist(solution$getValue(w))

#  Weighted random forest predictions
rf_predictions_weighted <- as.data.frame(predict(
    rf_model,
    test_data,
    predict.all = TRUE
))

#  Multiply predictions by weights per row
rf_predictions_weighted <- as.matrix(rf_predictions_weighted) %*%
    as.numeric(solution$getValue(w))

# Calculate weights
w <- Variable(ncol(rf_residuals))
objective <- quad_form(w, cov_matrix_shrinkage)
# objective <- norm2((t(w) %*% mean_vector)) + quad_form(w, cov_matrix_shrinkage)
constraints <- list(sum(w) == 1)
prob <- Problem(Minimize(objective), constraints)
solution <- solve(prob)
hist(solution$getValue(w))

#  Weighted random forest predictions
rf_predictions_weighted_shrinkage <- as.data.frame(predict(
    rf_model,
    test_data,
    predict.all = TRUE
))

#  Multiply predictions by weights per row
rf_predictions_weighted_shrinkage <- as.matrix(rf_predictions_weighted_shrinkage) %*% as.numeric(solution$getValue(w))

#  RMSE
rmse(mean(train_data$target), test_data$target) # naive benchmark
rmse(lm_predictions, test_data$target) # linear regression
rmse(rf_predictions[, 1], test_data$target) # random forest 1/n weights
rmse(rf_predictions_weighted, test_data$target) # random forest weighted
rmse(rf_predictions_weighted_shrinkage, test_data$target) # random forest weighted shrinkage
