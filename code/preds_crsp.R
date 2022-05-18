# load libraries
library(lubridate)
library(rlist)
library(doParallel)
library(dplyr)
library(glmnet)
library(xgboost)
library(ranger)
library(caret)
library(parallel)

#detect and set ncores 
ncores <- detectCores() - 1

# set working directory
setwd("~/Documents/Studium/PhD/Forecast_combination")

# settings
rolling_window <- 10 # in years
horizon <- 1
num_lags <- 12

# set up glmnet grid
grid_glmnet <- expand.grid(
  alpha = seq(0,1,0.1),
  lambda = seq(0.00001, 1, length = 50))

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

# read daily data
data = read.csv("data/CRSP/CRSPmonthly.csv", stringsAsFactors = FALSE)
data_factors = read.csv("data/CRSP/FFmonthly.csv", stringsAsFactors = FALSE)

# choose subset of columns and filter NAs
data = data[,c("PERMNO", "date", "RET")]
# length(unique(data$PERMNO)) # 24533 companies in total
data$date = round_date(as.Date(as.character(data$date), "%Y%m%d"),
                       unit = "month")
data = data[!(data$RET=="C" | as.character(data$RET)==""),]

# choose subset of columns and filter NAs
# data_factors$Date = as.Date(as.character(data_factors$Date), "%Y%m%d")
data_factors$Date = as.Date(paste0(data_factors$Date, "01"), "%Y%m%d") # for monthly data

# get unique comapny numbers
permno = unique(data$PERMNO)
permno <- permno

# # define function to get predictions per vintage
vintages <- sort(unique(data$date))
vintage_function <- function(vintage){
# for (vintage in as.list(vintages)) {
  vintage = as.Date(vintage)
  dates_bench = subset(data_factors$Date, data_factors$Date>=vintage & 
                         data_factors$Date<=vintage+years(rolling_window))
  # set up predictions dataframe
  predictions_vintage <- data.frame(matrix(ncol = nrow(grid_glmnet) +
                                             nrow(grid_rf)+nrow(grid_xgboost)+2
                                           ,nrow = 0))
  x <- c("Date", "permono", 
         paste0("glmnet_", seq(1:nrow(grid_glmnet))),
         paste0("xgboost_", seq(1:nrow(grid_xgboost))),
         paste0("rf_", seq(1:nrow(grid_rf))))
  colnames(predictions_vintage) <- x
  
  # set iteraator for prediciton saving
  k=1
  for (j in permno) {
    # print permno
    print(j)
    
    # extract data accodring to permno
    dates = data$date[data$PERMNO==j]
    dates_subset = subset(dates, dates>=vintage & dates<=as.Date(vintage)+years(rolling_window))
    
    if(length(dates_subset)==length(dates_bench)){
      # set iterator for dataframe
      
      # get returns and factors according to vintage
      returns = data[data$PERMNO==j,]
      returns = as.numeric(returns$RET[returns$date %in% dates_subset])
      factors = data_factors[data_factors$Date %in% dates_subset, 2:6]/100
      
      # get lagged features
      returns_lagged <- cbind(sapply(0:(num_lags-1), lag, x=returns))
      factors_lagged <- cbind(as.data.frame(lapply(0:(num_lags-1), lag, x=factors)))
      
      # bind all the features
      x <- cbind(returns_lagged, factors_lagged)

      # define target
      y <- returns
      
      # bind x and y and get training and oos ds
      ds <- cbind(y,x)
      training_ds <- ds[-nrow(ds), ]
      y_true <- ds[nrow(ds), "y"]
      
      # lead target for forecasting exercise
      training_ds$y <- lead(training_ds$y, horizon)
      
      # extract last observation for final forecast
      x_last_obs <- select(training_ds, -"y")
      x_last_obs <- x_last_obs[nrow(x_last_obs), ]
      
      # filter complete cases
      training_ds <- training_ds[complete.cases(training_ds),]
      
      # set up empty glm vector
      preds_glm <- vector()
      
      # produce forecasts
      for (i in 1:nrow(grid_glmnet)) {
        glm_fit <- glmnet(x = as.matrix(select(training_ds, -"y")), 
                          y = select(training_ds, "y")$y, 
                          alpha = grid_glmnet$alpha[i],
                          lambda = grid_glmnet$lambda[i])
        preds_glm[i] <- predict(glm_fit, as.matrix(x_last_obs))
      }

      # set up empty xgb vector
      preds_xgb <- vector()

      # produce forecasts
      for (i in 1:nrow(grid_xgboost)) {
        xgb_fit <- xgboost(data = as.matrix(select(training_ds, -"y")),
                           label = select(training_ds, "y")$y,
                           nrounds = grid_xgboost$nrounds[i],
                           eta = grid_xgboost$eta[i],
                           max_depth = grid_xgboost$max_depth[i],
                           gamma = grid_xgboost$gamma[i],
                           colsample_bytree = grid_xgboost$colsample_bytree[i],
                           min_child_weight = grid_xgboost$min_child_weight[i],
                           subsample = grid_xgboost$subsample[i],
                           verbose = 0)
                           # nthread = ncores)
        preds_xgb[i] <- predict(xgb_fit, as.matrix(x_last_obs))
      }

      # set up empty rf vector
      preds_rf <- vector()

      # produce forecasts
      for (i in 1:nrow(grid_rf)) {
        rf_fit <- ranger(x = as.matrix(select(training_ds, -"y")),
                         y = select(training_ds, "y")$y,
                         num.trees = grid_rf$num.trees[i],
                         mtry = grid_rf$mtry[i])
                         # num.threads = ncores)
        preds_rf[i] <- predict(rf_fit, as.matrix(x_last_obs))$predictions
      }
      predictions_vintage[k, ] <- c(paste0(max(as.Date(dates_subset))), 
                                    permno[j], preds_glm, preds_xgb, preds_rf)
      k=k+1
    }
  }
  print(vintage)
  save(predictions_vintage, file=paste0("data/CRSP/vintages_monthly/", 
                                        max(dates_bench),".RData"))
}

# # run function in parallel with mclapply
# # vintages = as.character(seq(from = as.Date("1964-01-01"), to = as.Date("2015-01-01"), by= "year"))
# vintages <- unique(data$date)[unique(data$date) >= as.Date("1980-01-01")]
mclapply(vintages, vintage_function, mc.cores = 1)
