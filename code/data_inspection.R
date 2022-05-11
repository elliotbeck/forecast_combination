# set wd
setwd("~/Documents/Studium/PhD/Forecast_combination")

# load library
library(dplyr)

# set startdate
startdate <- as.Date("1980-01-01")

# set wd
setwd("~/Documents/Studium/PhD/Forecast_combination")

# get vintages
vintages <- list.files(path="data/vintages")

# load names
names <- get(load("data/intersect_names.RData"))

# create list with predictor names
pred_names <- list()

# set up iterator 
i = 1

# loop through vintages
for(vintage in vintages){
  # load the data
  path_vintage <- paste0("data/vintages/", vintage)
  rawdata <- read.csv2(path_vintage, sep = ",")
  tcode <- rawdata[1,]
  
  # construct rawdata
  rawdata <- rawdata[-1,]
  
  # filter rows < 1980-01-01
  rawdata$sasdate <- as.Date(rawdata$sasdate, format = "%m/%d/%Y")
  rawdata <- rawdata[rawdata$sasdate > startdate, ]
  
  # make numeric
  rawdata[, 2:ncol(rawdata)] <- apply(rawdata[, 2:ncol(rawdata)], 
                                      as.numeric, MARGIN = 2)
  
  
  # get rid of columns/rows containing NAs
  rawdata <- rawdata[rowSums(is.na(rawdata))!=ncol(rawdata), ]
  rawdata <- rawdata[, colSums(is.na(rawdata))==0]
  
  # save in list
  pred_names[[i]] <- colnames(rawdata)
  
  # update iterator 
  i = i+1
}


intersect_names <- Reduce(intersect, pred_names)

# save the names
save(intersect_names, file="data/intersect_names.RData")

