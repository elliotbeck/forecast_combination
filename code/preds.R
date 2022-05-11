# set wd
setwd("~/Documents/Studium/PhD/Forecast_combination")

# load library
library(dplyr)
library(xgboost)
library(ranger)
library(stringr)
library(xts)
library(tsbox)
library(caret)
source("code/fredmd.R")

# settings
horizon <- 1
target <- "CPIAUCSL"

# set up glmnet grid
grid_glmnet <- expand.grid(
  alpha = seq(0,1,0.1),
  lambda = seq(0.00001, 2, length = 100))

# set up xgboost grid
grid_xgboost <- expand.grid(
  nrounds = seq(from = 100, to = 500, by = 100),
  eta = c(0.025, 0.1, 0.3),
  max_depth = seq(from = 2, to = 20, by = 4),
  gamma = 0,
  colsample_bytree = seq(from = 0.5, to = 0.9, by = 0.2),
  min_child_weight = 1,
  subsample = 1
)

# my tuneGrid object:
grid_rf <- expand.grid(
  num.trees = c(100, 500),
  mtry = seq(1,21,2),
  min.node.size = seq(1,6,1)
)

# set startdate
startdate <- as.Date("1980-01-01")

# get vintages
vintages <- list.files(path="data/vintages")

# load names
names <- get(load("data/intersect_names.RData"))

# create empty predictions vector
predictions <- data.frame(matrix(NA, nrow = length(vintages), 
                                 ncol = nrow(grid_glmnet)+nrow(grid_rf)+nrow(grid_xgboost)+1))

# set up iterator 
j = 1

# loop through vintages
for(vintage in vintages){
  # print iteration 
  print(vintage)
  
  # load the data
  path_vintage <- paste0("data/vintages/", vintage)
  rawdata <- read.csv2(path_vintage, sep = ",")
  
  # filter the columns that remain in whole period
  rawdata <- rawdata[, names]
  
  # extract transformation codes and remove first row
  tcode <- rawdata[1, 2:ncol(rawdata)]
  rawdata <- rawdata[-1,]
  
  # change tcode of target to mom changes (case 8)
  tcode[names(tcode)==target] <- 8
  
  # filter rows < startdate
  rawdata$sasdate <- as.Date(rawdata$sasdate, format = "%m/%d/%Y")
  rawdata <- rawdata[rawdata$sasdate > as.Date("1980-01-01"), ]
  
  # make numeric
  rawdata[, 2:ncol(rawdata)] <- apply(rawdata[, 2:ncol(rawdata)], 
                                      as.numeric, MARGIN = 2)
  
  
  # get rid of rows containing NAs
  rawdata <- rawdata[rowSums(is.na(rawdata))!=ncol(rawdata), ]

  # transform data
  data <- fredmd(rawdata, tcode)
  
  # remove rows with NAs due to transformations
  data <- data[complete.cases(data), ]
  
  # get target
  y <- select(data, all_of(target))
  
  # lead target for forecasting exercise
  y <- na.omit(lead(y, horizon))
  
  # get predictors
  predictors <- select(data, -c("sasdate"))
  
  # get last observation and predictors for training (x)
  x <- predictors[1:nrow(y), ]
  x_lastobs <- predictors[nrow(y)+horizon, ]
  
  # set up empty glm vector
  preds_glm <- vector()
  
  # produce forecasts
  for (i in 1:nrow(grid_glmnet)) {
    glm_fit <- glmnet(x, y$CPIAUCSL, 
                      alpha = grid_glmnet$alpha[i],
                      lambda = grid_glmnet$lambda[i])
    preds_glm[i] <- predict(glm_fit, as.matrix(x_lastobs))
  }
  
  # set up empty xgb vector
  preds_xgb <- vector()

  # produce forecasts
  for (i in 1:nrow(grid_xgboost)) {
    xgb_fit <- xgboost(as.matrix(x), y$CPIAUCSL,
                       nrounds = grid_xgboost$nrounds[i],
                       eta = grid_xgboost$eta[i],
                       max_depth = grid_xgboost$max_depth[i],
                       gamma = grid_xgboost$gamma[i],
                       colsample_bytree = grid_xgboost$colsample_bytree[i],
                       min_child_weight = grid_xgboost$min_child_weight[i],
                       subsample = grid_xgboost$subsample[i], 
                       verbose = 0)
    preds_xgb[i] <- predict(xgb_fit, as.matrix(x_lastobs))
  }
  
  # set up empty glm vector
  preds_rf <- vector()
  
  # produce forecasts
  for (i in 1:nrow(grid_rf)) {
    rf_fit <- ranger(x = as.matrix(x), y = y$CPIAUCSL, 
                     num.trees = grid_rf$num.trees[i],
                     mtry = grid_rf$mtry[i])
    preds_rf[i] <- predict(rf_fit, as.matrix(x_lastobs))$predictions
  }

  # save prediction
  predictions[j,] <- c(paste0(str_remove(vintage, ".csv"), "-01"), 
                       preds_glm, preds_rf, preds_xgb)
  
  # update iterator 
  j = j+1
}

# save predictions 
save(predictions, file="data/predictions.RData")

