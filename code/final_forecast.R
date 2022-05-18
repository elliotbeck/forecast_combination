# set wd
setwd("~/Documents/Studium/PhD/Forecast_combination")

# load libraries
library(lubridate)
library(xts)
library(tsbox)
source("code/fredmd.R") 
source("code/rmse.R")
source("code/qis.R")

# settings
horizon <- 1
target <- "CPIAUCSL"
names <- get(load("data/intersect_names.RData"))

# load results and cpi
load("data/predictions.RData")
rawdata <- read.csv("data/vintages/2022-04.csv")

# load weights
load("~/Documents/Studium/PhD/Forecast_combination/data/weights_norm1_1.5.RData")

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

# assure predictions are numeric
predictions[, 2:ncol(predictions)] <- apply(predictions[, 2:ncol(predictions)], 
                                            as.numeric, MARGIN = 2)

# delete last prediction, not yet observed...
predictions <- predictions[-nrow(predictions), ]
predictions_xts <- xts(predictions[, 2:ncol(predictions)],
                       order.by = as.Date(predictions$X1))

# bring to xts format
cpi_xts <- window(xts(data$CPIAUCSL, order.by = as.Date(data$sasdate)), 
                  start = index(predictions_xts)[1],
                  end = index(predictions_xts)[nrow(predictions_xts)])


# sanity check
weights <- weights[, colSums(is.na(weights))!=nrow(weights)]

# create final forecast xts
final_forecast <- xts(rep(NA, ncol(weights)), 
                      order.by = as.Date(colnames(weights))+months(1))

# get final forecasts
for (i in 1:(ncol(weights)-1)) {
  date <- as.Date(colnames(weights)[i]) + months(1)
  preds <- predictions_xts[index(predictions_xts) == date, ]
  final_forecast[i] <- sum(weights[,i] * preds, na.rm=TRUE)
}
final_forecast <- na.omit(final_forecast)

# create mean forecast xts
mean_forecast <- xts(rep(NA, ncol(weights)), 
                      order.by = as.Date(colnames(weights))+months(1))

# get mean forecasts
for (i in 1:(ncol(weights)-1)) {
  date <- as.Date(colnames(weights)[i]) + months(1)
  preds <- predictions_xts[index(predictions_xts) == date, ]
  n <- sum(as.numeric(!is.na(weights[,i])))
  mean_forecast[i] <- sum(as.numeric(!is.na(weights[,i])) * preds, na.rm=TRUE)/n
}
mean_forecast <- na.omit(mean_forecast)

# naive mean
dates <- index(cpi_xts[index(cpi_xts)>=index(final_forecast)[1] &
                   index(cpi_xts)<=index(final_forecast)[nrow(final_forecast)]])
mean_forecast_naive <- xts(rep(NA, length(dates)), order.by = as.Date(dates))
i=1
for (date in as.list(dates)) {
  preds <- predictions_xts[index(predictions_xts) == date, ]
  mean_forecast_naive[i] <- mean(as.numeric(preds))
  i=i+1
}

# get rmses
cpi_xts <- cpi_xts[index(cpi_xts)>=index(final_forecast)[1] &
                     index(cpi_xts)<=index(final_forecast)[nrow(final_forecast)]]

# include drift term in rw, as cpi series has non zero mean
rw <- cumsum(cpi_xts)/cumsum(rep(1,nrow(cpi_xts)))

ts_plot(cbind(cpi_xts, final_forecast, mean_forecast, 
              mean_forecast_naive, rw))

f_rmse(rw, cpi_xts)
f_rmse(final_forecast, cpi_xts)
f_rmse(mean_forecast, cpi_xts)
f_rmse(mean_forecast_naive, cpi_xts)

