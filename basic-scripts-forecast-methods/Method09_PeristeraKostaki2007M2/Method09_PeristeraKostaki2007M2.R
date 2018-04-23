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
#    Method08_PeristeraKostaki2007M2.R
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
#   1. Peristera P, Kostaki A (2007) Modeling fertility in modern populations. Demographic Research 16(6):141-194.
#       
#       
###############################################################################
library(forecast)
##############################################################################
#Peristera and Kostaki (Model 2):
################################ 
PK2_fx <- function(par,x){
  c1 <- par[1]
  c2 <- par[2]
  mu1 <- par[3]
  mu2 <- par[4]
  d11 <- par[5]
  d12 <- par[6]
  d2 <- par[7]
  fx1 <- c1*exp(-((x[which(x<=mu1)]-mu1)/d11)^2)
  fx2 <- c1*exp(-((x[which(x>mu1)]-mu1)/d12)^2)
  fx3 <- c2*exp(-((x-mu2)/d2)^2)
  fx <- c(fx1,fx2)+fx3
  names(fx) <- x
  return(fx)
}

Fit_PK2_fx <- function(par,data,method,x){
  fctn_PK2 <- function(par, data){
    res <- sum((data - PK2_fx(par,x))^2)
    return(res)
  }
  fit <- optim(par,fctn_PK2,gr=NULL,data,method=method,control=list(trace=F)) 
  return(fit$par)
}


if(0>1){
  aut.complete = ASFR[,"2014"]
  fitPK2.complete <- Fit_PK2_fx(par=c(0.1,0.1,20,28,3,10,7),data=aut.complete,x=15:44,method="Nelder-Mead") 
  fitPK2.complete 
  
  plot(x=15:44,y=aut.complete,col="red")
  lines(x=15:44,y=PK2_fx(fitPK2.complete,15:44),col="red")
  
}

###########################################
Method09_PeristeraKostaki2007M2.R <- function(ASFR,
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
                method="Peristera Kostaki M2 (2007)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("PeristeraKostakiM2",parameter,sep="_"),"PeristeraKostakiM2"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  
  if(parameter=="arima" | parameter=="freezed"){
    ######## get paramerters
    param <- matrix(NA,nrow=7,ncol=(year3-year1+1))
    colnames(param) <- paste(year1 + c(1:ncol(param))-1)
    param[,colnames(fert.raw)] <- apply(fert.raw,2,Fit_PK2_fx,par=c(0.1,0.1,20,28,3,10,7),method="Nelder-Mead",x=age1:age2)
    
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
      param[4,] <- projection.of.knots.arima(param[4,])
      param[5,] <- projection.of.knots.arima(param[5,])
      param[6,] <- projection.of.knots.arima(param[6,])
      param[7,] <- projection.of.knots.arima(param[7,])
    }
    if(parameter[1]=="freezed"){
      param[1,] <- projection.of.knots.freezed( param[1,])
      param[2,] <- projection.of.knots.freezed( param[2,])
      param[3,] <- projection.of.knots.freezed( param[3,])
      param[4,] <- projection.of.knots.freezed( param[4,])
      param[5,] <- projection.of.knots.freezed( param[5,])
      param[6,] <- projection.of.knots.freezed( param[6,])
      param[7,] <- projection.of.knots.freezed( param[7,])
    }
    
    ##### project of profile
    for(yr in colnames(raw.data[,is.na(raw.data[1,])])){
      raw.data[,yr] <- PK2_fx(param[,yr],age1:age2)
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
                method="Peristera Kostaki M2 (2007)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("PeristeraKostakiM2",parameter,sep="_"),"PeristeraKostakiM2"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  }
  
  ######## The original method of PK2
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
        
        #PK2:
        pk2 <-  Fit_PK2_fx(par=c(0.1,0.1,20,28,3,10,7),data=data,x=x,method="Nelder-Mead") 
        pk2.est <- PK2_fx(pk2,age1:age2)
        est.cdata[,coh] <- pk2.est	
        
      }
    }
    
    ######## original
    ########
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
                predASFR=asfr_cohort_to_period(est.cdata),
                predCASFR=est.cdata,
                method="Peristera Kostaki M2 (2007)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("PeristeraKostakiM2",parameter,sep="_"),"PeristeraKostakiM2"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  }
  
}
