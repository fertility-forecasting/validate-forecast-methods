#######################################################################
#######################################################################
## 'validate-forecast-methods' is a program that provides the R source code
## used for the calculations in the project: 
## Bohk-Ewald, Christina, Peng Li, and Mikko Myrskylä (2017). 
## Assessing the accuracy of cohort fertility forecasts. 
## Presented in session: Statistical methods in demography 
## at the PAA 2017 Annual Meeting, Chicago, IL, USA, April 27-April 29, 2017. 
## (c) Copyright 2018, Christina Bohk-Ewald, Peng Li, Mikko Myrskylä

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program (see LICENSE.txt in the root of this source code).
## If not, see <http://www.gnu.org/licenses/>.
#######################################################################
#######################################################################

###############################################################################
#    Project: Cohort fertility forecastes and their accuracy
#    CBE, PL and MM
#    2017.05.19
#    
#    Method:   
#    Method03_CoaleTrussell1974.R
############################################################
#### Input file: 
#     ASFR: age-specific fertility rate
#     
#     
#### Parameters
#      joy: jump of year
#      obs: length of observation
#      age1 : star age
#      age2 : end age
#      parameter : parameter of the method
#                projection method of parameters (alpha, beta) [arima/freezed]
#                arima: max ordre of (p,d,q) is (5,2,2) using auto.arima
#                freezed: use the last observed parameter (alpha, beta) 
#      len : length of forecasting period
#      pop : population
#### Output: 
#      obsASFR: observed ASFR, age * period
#      obsCASFR: observed ASFR, age * cohort
#      predASFR: forecasted ASFR, age * period
#      predCASFR: forecasted ASFR, age * cohort
#      method: name of the method
#      parameter: parameter of the method
#      year: starting year of observed period, JOY and last year of forecasting
#      age: range of age
#      obs: number of years of observation
#      len: forecasting length
#      cohort: c1: oldest cohort at the starting year of observation
#              c2: youngest cohort at the starting year of observation
#              c3: oldest cohort at JOY
#              c4: youngest cohort at JOY
#### NOTE:
#      
#### Ref:
#   1. Coale AJ, Trussell TJ (1974) Model Fertility Schedules: Variations in The Age Structure of
#      Childbearing in Human Populations. Population Index 40(2):185-258.
#
###############################################################################
library(forecast)
##############################################################################
#Coale Trussell 1974 method
################################ 
#### use 15:49
CT <- function (M = 1, m = 0,age1=15,age2=49) 
{
  nf <- c(0.325, 0.375, 0.421, 0.46, 0.475, 
          0.477, 0.475, 0.47, 0.465, 0.46, 0.455, 0.449, 0.442, 
          0.435, 0.428, 0.42, 0.41, 0.4, 0.389, 0.375, 0.36, 0.343, 
          0.325, 0.305, 0.28, 0.247, 0.207, 0.167, 0.126, 0.087, 
          0.055, 0.035, 0.021, 0.011, 0.003)
  rn <- c(0, 0, 0, 0, 0, 0.004, 0.03, 0.06, 0.1, 0.15, 
          0.2, 0.25, 0.31, 0.37, 0.44, 0.52, 0.6, 0.68, 0.76, 0.83, 
          0.9, 0.97, 1.04, 1.11, 1.18, 1.25, 1.32, 1.39, 1.46, 
          1.53, 1.59, 1.64, 1.67, 1.69, 1.7) * (-1)
  ages <- 15:49
  rval <- M * nf * exp(m * rn)
  names(rval) <- ages
  return(rval[paste(age1:age2)])
}
#### using optim to estimate parameters
fitCT_optim <- function (initialpar = c(1, 1), data, age1=15,age2=49, Method = "Nelder-Mead",...) 
{
  fCT <- function(z, rASMFR,age1,age2) {
    nc <- length(rASMFR)
    mASMFR <- CT(z[1], z[2],age1,age2)
    res <- sum((rASMFR - mASMFR)^2)
    RMSE <- sqrt(res/nc)
    return(RMSE)
  }
  rCT <- optim(initialpar, fCT, gr = NULL, data,age1,age2, method = Method, ...)
  return(c(rCT$par, rCT$value, rCT$convergence))
}
#res <- fitCT_optim(c(1,1),asfr,age1=15,age2=44)

#### using linear model to estimate parameters

fitCT_linear <- function (data,age1=15,age2=49) 
{
  nf <- c(0.325, 0.375, 0.421, 0.46, 0.475, 
          0.477, 0.475, 0.47, 0.465, 0.46, 0.455, 0.449, 0.442, 
          0.435, 0.428, 0.42, 0.41, 0.4, 0.389, 0.375, 0.36, 0.343, 
          0.325, 0.305, 0.28, 0.247, 0.207, 0.167, 0.126, 0.087, 
          0.055, 0.035, 0.021, 0.011, 0.003)
  rn <- c(0, 0, 0, 0, 0, 0.004, 0.03, 0.06, 0.1, 0.15, 
          0.2, 0.25, 0.31, 0.37, 0.44, 0.52, 0.6, 0.68, 0.76, 0.83, 
          0.9, 0.97, 1.04, 1.11, 1.18, 1.25, 1.32, 1.39, 1.46, 
          1.53, 1.59, 1.64, 1.67, 1.69, 1.7) * (-1)

  obj <- lm(log(data/nf[c(age1:age2)-15+1])~rn[c(age1:age2)-15+1])
  M <- exp(obj$coef[1])
  m <- obj$coef[2]
  return(c(M,m))
}

#### using 20-24 to estimate M, 25-49 to estimate m
#data <-  c(0, 0, 0, Jfert$ASMFR2000[1:35])
fitCT_twostep <- function (data,age1=15,age2=49) 
{
  nf <- c(0.325, 0.375, 0.421, 0.46, 0.475, 
          0.477, 0.475, 0.47, 0.465, 0.46, 0.455, 0.449, 0.442, 
          0.435, 0.428, 0.42, 0.41, 0.4, 0.389, 0.375, 0.36, 0.343, 
          0.325, 0.305, 0.28, 0.247, 0.207, 0.167, 0.126, 0.087, 
          0.055, 0.035, 0.021, 0.011, 0.003)
  rn <- c(0, 0, 0, 0, 0, 0.004, 0.03, 0.06, 0.1, 0.15, 
          0.2, 0.25, 0.31, 0.37, 0.44, 0.52, 0.6, 0.68, 0.76, 0.83, 
          0.9, 0.97, 1.04, 1.11, 1.18, 1.25, 1.32, 1.39, 1.46, 
          1.53, 1.59, 1.64, 1.67, 1.69, 1.7) * (-1)
  
  M <- sum(data[c(20:24)-age1+1])/sum(nf[c(20:24)-15+1]) 
  m <- mean(log(data[c(25:age2)-age1+1]/nf[c(25:age2)-15+1])/rn[c(25:age2)-15+1])
  
  return(c(M,m))
}
##############################################
fitCT_optim_pred <- function(asfr,age1=15,age2=49,method="optim"){
  if(method=="optim"){
    res <- fitCT_optim(c(1,1),asfr,age1,age2)
    FLAG <- res[4]
    while (FLAG>0) {
      res <- fitCT_optim(res[1:2], asfr)
      FLAG <- res[4]
    }
    M <- res[1]
    m <- res[2]
  }
  
  if(method=="linear"){
    obj <- fitCT_linear(asfr,age1,age2)
    M <- obj[1]
    m <- obj[2]
  }
  
  if(method=="twostep"){
    obj <- fitCT_twostep(asfr,age1,age2)
    M <- obj[1]
    m <- obj[2]
  }
  return(c(M,m))
  #print(c(M,m))
  #fitted <- CT(M,m)
  #return(fitted)
}

###########################################
Method03_CoaleTrussell1974.R <- function(ASFR,
                                         joy = 1985, 
                                         obs = 30,
                                         age1 = 15, 
                                         age2 = 49,
                                         parameter = c("linear","arima"),
                                         len = 50,
                                         pop=""
){
  
  ################ Data step
  if(is.na(joy))
    joy <- max(as.numeric(colnames(ASFR)[!is.na(apply(ASFR,2,sum))]))
  ########
  year1 <- joy - obs + 1
  year2 <- joy
  year3 <- year2 + len
  ########
  # extract input data
  fert.raw <- ASFR[paste(age1:age2),paste(year1:year2)]
  fert.raw[fert.raw<=0.00001] <- 0.00001
  raw.data <- matrix(NA,nrow=nrow(fert.raw),ncol=(year3-year1+1))
  colnames(raw.data) <- paste(year1 + c(1:ncol(raw.data))-1)
  rownames(raw.data) <- rownames(fert.raw)
  raw.data[,colnames(fert.raw)] <- fert.raw

  ########
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))]))))
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
                predASFR=NA,
                predCASFR=NA,
                method="Coale Trussell (1974)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("CoaleTrussell",paste(parameter,collapse = "_"),sep="_"),"CoaleTrussell"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  ######## get paramerters
  param <- matrix(NA,nrow=2,ncol=(year3-year1+1))
  colnames(param) <- paste(year1 + c(1:ncol(param))-1)
  param[,colnames(fert.raw)] <- apply(fert.raw,2,fitCT_optim_pred,age1,age2,parameter[1])
  
  ######## projections of knots 
  #### arima model
  projection.of.knots.arima <- function(x){
    obj <- auto.arima(as.ts(x[!is.na(x)]),max.p=5,max.d = 2,max.q = 2)
    pre <- forecast(obj,h = length(x)-max(which(!is.na(x))))
    x[c(max(which(!is.na(x)))+1):length(x)] <- pre$mean
    return(x)
  }
  #### freezed rate
  projection.of.knots.freezed <- function(x){
    ind <- !is.na(x)
    x[max(which(ind)+1):length(x)] <- x[max(which(ind))]
    return(x)
  }
  
  #### prejection of parameter
  if(parameter[2]=="arima"){
    param[1,] <- projection.of.knots.arima(param[1,])
    param[2,] <- projection.of.knots.arima(param[2,])
  }
  if(parameter[2]=="freeze"){
    param[1,] <- projection.of.knots.freezed( param[1,])
    param[2,] <- projection.of.knots.freezed( param[2,])
  }
  
  ##### project of profile
  for(yr in colnames(raw.data[,is.na(raw.data[1,])])){
    raw.data[,yr] <- CT(param[1,yr],param[2,yr],age1,age2)
  }
  
  ################ Output step
  ####predASFR=raw.data[,paste(year2+1:len)] ##this only output the predicted part
  predASFR=raw.data                          ##this output observed and predicted part
  predCASFR <- asfr_period_to_cohort(raw.data)
  #predCV <- predCASFR[,paste((year2-age1):(year2-age2+1))]
  
  ########
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              method="Coale Trussell (1974)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("CoaleTrussell",paste(parameter,collapse = "_"),sep="_"),"CoaleTrussell"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}
