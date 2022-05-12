# set wd
setwd("~/Documents/Studium/PhD/Forecast_combination")

# load libraries
library(CVXR)
library(HDShOP)
library(lubridate)
library(xts)
source("code/fredmd.R") 
source("code/qis.R")

# settings
horizon <- 1
target <- "CPIAUCSL"
names <- get(load("data/intersect_names.RData"))
rolling_window <- 5

# load results and cpi
load("data/predictions.RData")
rawdata <- read.csv("data/vintages/2022-04.csv")
load("data/weights_final.RData")

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

# create sequence of dates
vintages <- seq(as.Date("2010/1/1"), as.Date("2022/3/1"), "months")

# create empty weights dataframe
weights <- as.data.frame(matrix(rep(NA, ncol(predictions_xts)*length(vintages)), 
                                ncol = length(vintages)))
rownames(weights) <- colnames(predictions_xts)
colnames(weights) <- as.character(vintages)
j=1
for (vintage in as.list(vintages)) {
  # get status updates
  print(vintage)
  
  # create rolling window
  preds_vintage <- predictions_xts[index(predictions_xts)<=vintage & 
                             index(predictions_xts)>=vintage-years(rolling_window), ]
  cpi_vintage <- cpi_xts[index(cpi_xts)<=vintage & 
                            index(cpi_xts)>=vintage-years(rolling_window), ]
  
  # calculate errors
  errors_vintage <- preds_vintage - as.numeric(cpi_vintage)
  
  # kick out super correlated predictions (>95%)
  # cor_errors_temp <- cov2cor(qis(errors))
  cor_errors_temp <- cor(errors_vintage)
  cor_errors_temp[upper.tri(cor_errors_temp)] <- 0
  diag(cor_errors_temp) <- 0
  while(max(cor_errors_temp)>0.95){
    # get maximum name
    max_cor <- which(cor_errors_temp == max(cor_errors_temp), arr.ind = TRUE)

    # calculate without maximum correlation column
    errors_vintage <- errors_vintage[, -max_cor[2]]
    # cor_errors_temp <- cov2cor(qis(errors))
    cor_errors_temp <- cor(errors_vintage)
    cor_errors_temp[upper.tri(cor_errors_temp)] <- 0
    diag(cor_errors_temp) <- 0
  }

  # get names of selected predictions
  names_selected <- colnames(cor_errors_temp)
  length(names_selected)
  # names_selected <- rownames(weights_final)[!is.na(weights_final[, colnames(weights_final)==vintage])]
  
  # extract selected erros and calculate covariance and sample mean
  preds_selected <- preds_vintage[, names_selected]
  errors_selected <- preds_selected - as.numeric(cpi_vintage)
  # cov_selected <- qis(errors_selected) # ledoit wolf non linear shrink
  cov_selected <- cov(errors_selected) # sample cov
  # mean_selected <- as.matrix(mean_bs(t(errors_selected))$mean) # bayes stein
  mean_selected <- colMeans(errors_selected) # sample mean
  
  # convex optimization
  w <- Variable(ncol(errors_selected))
  objective <- norm2((t(w)%*%mean_selected)) + quad_form(w, cov_selected)
  # constraints <- list(sum(w) == 1, w>=0)
  constraints <- list(sum(w) == 1, w>=0)
  prob <- Problem(Minimize(objective), constraints)
  solution <- solve(prob)
  
  # add calculated weights to df
  weights[rownames(weights)%in%names_selected, j] <- solution$getValue(w)
  
  j=j+1
}

# save predictions 
save(weights, file="data/weights_rolling5.RData")
