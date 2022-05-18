# load libraries
library(lubridate)
library(rlist)
library(sandwich)
library(doParallel)


# set working directory
setwd("~/Documents/Studium/PhD/Forecast_combination")

# read daily data
data = read.csv("~/Desktop/CRSPdaily.csv", stringsAsFactors = FALSE)
data_factors = read.csv("~/Desktop/FFdaily.csv", stringsAsFactors = FALSE)

# choose subset of columns and filter NAs
data = data[,c("PERMNO", "date", "RET")]
length(unique(data$PERMNO)) # 24533 companies in total
data$date = as.Date(as.character(data$date), "%Y%m%d")
data = data[!(data$RET=="C" | as.character(data$RET)==""),]

# choose subset of columns and filter NAs
data_factors$Date = as.Date(as.character(data_factors$Date), "%Y%m%d")

# load als function
load("WLS/wls.als.RData")

# define function to estimate ols regressions
permno = unique(data$PERMNO)
vintage_function <- function(vintage){
  vintage = as.Date(vintage)
  dates_bench = subset(data_factors$Date, data_factors$Date>=vintage & data_factors$Date<=vintage+years(5))
  regressions_res_factor <- list()
  
  for (j in permno) {
    dates = data$date[data$PERMNO==j]
    dates_subset = subset(dates, dates>=vintage & dates<=as.Date(vintage)+years(5))
    
    if(length(dates_subset)==length(dates_bench)){
      returns = data[data$PERMNO==j,]
      returns = as.numeric(returns$RET[returns$date %in% dates_subset])
      exc_returns = returns - data_factors$RF[data_factors$Date %in% dates_subset]
      factors = data_factors[data_factors$Date %in% dates_subset, 2:6]/100
      
      #get OLS HC results
      fit = lm(exc_returns~factors$SMB+factors$HML+factors$RMW+factors$CMA, x=TRUE, y=TRUE)
      sumry = as.data.frame(as.data.frame(summary(fit)[["coefficients"]])[,1])
      colnames(sumry) <- "ols_coef"
      sumry$ols_HC3_std = as.numeric(diag(vcovHC(fit)))
      
      # get ALS/WLS results
      result_als_wls = lm.wls.als(fit.ols=fit, spec.wool = TRUE)
      
      sumry$wls_coef <- summary(result_als_wls$fit.wls)$coefficients[,1]
      sumry$wls_HC3_std <- as.numeric(diag(vcovHC(result_als_wls$fit.wls)))
      
      sumry$als_coef <- summary(result_als_wls$fit.als)$coefficients[,1]
      sumry$als_HC3_std <- as.numeric(diag(vcovHC(result_als_wls$fit.als)))
      # # append results to list
      regressions_res_factor = list.append(regressions_res_factor, cbind(sumry, j))
      print(regressions_res_factor)
    }
  }
  print(vintage)
  save(regressions_res_factor, file=paste0("~/Desktop/R Code/Results/results_factors_", vintage,".RData"))
}

# run function in parallel with mclapply
vintages = as.character(seq(from = as.Date("1964-01-01"), to = as.Date("2015-01-01"), by= "year"))

mclapply(vintages, vintage_function, mc.cores = 1)