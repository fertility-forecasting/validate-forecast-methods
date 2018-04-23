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
#    Method02_CoaleMcNeil1972.R
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
#                projection method of parameters (alpha, beta) [original/arima/freezed]
#                original: for each cohort, use observed values to estimate parameters, then use this parameter to project
#                arima: max ordre of (p,d,q) is (5,2,2) using auto.arima
#                freeze: use the last observed parameter (alpha, beta) 
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
#   1. Coale A, McNeil D (1972) The distribution by age at ???rst marriage in a female cohort. Journal
#      of the American Statistical Association 67(340):743-749.
#      
#       
###############################################################################
#install.packages("nuptfer")
#library(devtools)
#install_github("grodri/nuptfer")
#install.packages("fmsb")
require(nuptfer)
require(fmsb)
require(forecast)
##############################################################################
#Coale Trussell 1974 method
################################ 
CoaleMcNeil_fx <- function (C = 0.8, a0 = 15, K = 2, ages=15:44){
  
  gs <- function(x) {
    0.1946 * exp(-0.174 * (x - 6.06) - exp(-0.2881 * (x -6.06)))
  }
  GCem <- GCpm <- GCfm <- C * sapply((ages - a0)/K, gs)
  for (i in 2:length(ages)) {
    GCpm[i] <- GCfm[i] * (1 - GCem[i - 1])
    GCem[i] <- GCpm[i] + GCem[i - 1]
  }
  mu <- sum(ages * GCpm)
  V <- sum((ages - mu)^2 * GCpm, na.rm = T)

  return(list(g = GCfm, G = GCem, mu = mu))
  #return(list(g = GCfm, G = GCem, mu = mu, sigma = sqrt(V)))
}

Fit_CoaleMcNeil_fx <- function (par = c(0.8, 18, 2), data, ages = 15:44, mode = 1,method = "Nelder-Mead"){
  if (mode == 1) {
    fCM <- function(z, Dz = data, Mz = mode, Az = ages) {
      nc <- length(Dz)
      Azz <- Az[1]:Az[nc] - 14
      mmz <- CoaleMcNeil_fx(z[1], z[2], z[3])$g[Azz]
      res <- sum((Dz - mmz)^2, na.rm = T)
      RMSE <- sqrt(res/nc)
      return(RMSE)
    }
  }
  else {
    fCM <- function(z, Dz = data, Mz = mode, Az = ages) {
      nc <- length(Dz)
      Azz <- Az[1]:Az[nc] - 14
      mmz <- CoaleMcNeil_fx(z[1], z[2], z[3])$G[Azz]
      res <- sum((Dz - mmz)^2, na.rm = T)
      RMSE <- sqrt(res/nc)
      return(RMSE)
    }
  }
  rCM <- optim(par, fCM, Dz = data, Mz = mode, Az = ages, method = method)
  return(c(rCM$par))
}

if(0>1){
  aut.complete = ASFR[,"2014"]
  res2 <- Fit_CoaleMcNeil_fx(par=c(1, 19, 1), data=aut.complete[as.character(15:44)], ages=15:44, mode=1, method="Nelder-Mead")
  fitted2 <- CoaleMcNeil_fx(res2[1], res2[2], res2[3],ages=15:44)
  plot(x=15:34,y=aut.complete[as.character(15:34)],xlim=c(15,44))
  lines(x=15:44,y=fitted2$g,col="blue")
  lines(x=15:44,y=aut.complete,col="red")
}

###########################################
Method02_CoaleMcNeil1972.R <- function(ASFR,
                                         joy = 1985, 
                                         obs = 30,
                                         age1 = 15, 
                                         age2 = 44,
                                         parameter = c("original"),
                                         len = 30,
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
                method="Coale McNeil (1972)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("CoaleMcNeil",parameter,sep="_"),"CoaleMcNeil"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  if(parameter=="arima" | parameter=="freezed"){
  ######## get paramerters
  param <- matrix(NA,nrow=3,ncol=(year3-year1+1))
  colnames(param) <- paste(year1 + c(1:ncol(param))-1)
  param[,colnames(fert.raw)] <- apply(fert.raw,2,Fit_CoaleMcNeil_fx,par=c(1,18,1),method="Nelder-Mead",ages=age1:age2)
  
  ######## projections of knots 
  #### arima model
  projection.of.knots.arima <- function(x){
    obj <- auto.arima(as.ts(x[!is.na(x)]),max.p=5,max.d = 2,max.q = 2)
    pre <- forecast(obj,h = length(x)-max(which(!is.na(x))))
    x[c(max(which(!is.na(x)))+1):length(x)] <- pre$mean
    x[x < 0.01] <- 0.01
    return(x)
  }
  #### freezed rate
  projection.of.knots.freezed <- function(x){
    ind <- !is.na(x)
    x[max(which(ind)+1):length(x)] <- x[max(which(ind))]
    return(x)
  }
  
  #### prejection of parameter
  if(parameter[1]=="arima"){
    param[1,] <- projection.of.knots.arima(param[1,])
    param[2,] <- projection.of.knots.arima(param[2,])
    param[3,] <- projection.of.knots.arima(param[3,])
  }
  if(parameter[1]=="freeze"){
    param[1,] <- projection.of.knots.freezed( param[1,])
    param[2,] <- projection.of.knots.freezed( param[2,])
    param[3,] <- projection.of.knots.freezed( param[3,])
  }
  
  ##### project of profile
  for(yr in colnames(raw.data[,is.na(raw.data[1,])])){
    raw.data[,yr] <- CoaleMcNeil_fx(param[1,yr],param[2,yr],param[3,yr],age1:age2)$g
  }
  
  ################ Output step
  predASFR=raw.data[,paste(year2+1:len)]
  predCASFR <- asfr_period_to_cohort(raw.data)
  #predCV <- predCASFR[,paste((year2-age1):(year2-age2+1))]
  ########
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              method="Coale McNeil (1972)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("CoaleMcNeil",parameter,sep="_"),"CoaleMcNeil"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  }
  
  ######## The original method of PK
  if(parameter == "original"){
    # age*period to age*cohort 
    raw.cdata <- asfr_period_to_cohort(fert.raw)
    est.cdata <- raw.cdata
    for(coh in colnames(raw.cdata)){
      data <- raw.cdata[,coh]
      ind <- !is.na(data)
      if(sum(ind)){
        data <- data[ind]
        x <- c(age1:age2)[ind]

        #CoaleMcNeil:
        cmn <-  Fit_CoaleMcNeil_fx(par=c(0.8, 20, 2),data=data,ages=x,mode=1,method="Nelder-Mead") 
        cmn.est <- CoaleMcNeil_fx(cmn[1],cmn[2],cmn[3],ages=age1:age2)
        est.cdata[,coh] <- cmn.est$g	
        
        
      }
    }
    
    ######## original
    ########
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
                predASFR=asfr_cohort_to_period(est.cdata),
                predCASFR=est.cdata,
                method="Coale McNeil (1972)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("CoaleMcNeil",parameter,sep="_"),"CoaleMcNeil"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  }

}
