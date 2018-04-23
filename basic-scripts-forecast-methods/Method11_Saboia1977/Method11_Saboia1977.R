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
#    2016.10.04
#    
#    Method: Saboia 1977  
#    Method11_Saboia1977.R
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
#      parameter : order
#             order : provide fixed numeric order; otherwise use auto.arima to select order automatically.
#             for auto.arima, max.p = 5,max.d=2,max.q=2
#      len : length of forecasting period
#      pop : population

#### Output: 
#      obsASFR: observed ASFR, age * period
#      obsCASFR: observed ASFR, age * cohort
#      predASFR: forecasted ASFR, age * period
#      predASFRlowerPI[80,90] : lower boundary of PI
#      predASFRupperPI[80/90] : upper boundary of PI
#      predCASFR: forecasted ASFR, age * cohort
#      predCASFRlowerPI[80,90] : lower boundary of PI
#      predCASFRupperPI[80/90] : upper boundary of PI
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
#
#### NOTE:
#
#### Ref:
#   1. Saboia JLM (1977) Autoregressive Integrated Moving Average (ARIMA) Models for Birth Forecasting. 
#      Journal of the American Statistical Association 72(358):264-270.     
#        
###############################################################################
#install.packages("forecast")
library(forecast)
############################################
#### ARIM model
tfr_arima_pred <- function(x,order = c(4,1,1),len=30,auto.bound=T){
  if(!is.numeric(order)){
    if(auto.bound)
      obj <- auto.arima(x,max.p=5,max.d=2,max.q=2)
    if(!auto.bound)
      obj <- auto.arima(x)
  }
  if(is.numeric(order) & sum(order>=0)==3){
    obj <- Arima(x, order = order)
  }
  
  #### prediction
  pred <- forecast(obj, len)
  
  #### outpute
  out <- rbind(pred$mean,
               pred$mean,
               pred$lower[,1],
               pred$upper[,1],
               pred$lower[,2],
               pred$upper[,2])
  row.names(out) <- c("TFR_ObsPre","TFR_FitPre","lowerPI80","upperPI80","lowerPI95","upperPI95")
  out <- cbind(matrix(rep(x,6),nrow=6,byrow = T),out)
  out[2,c(1:length(x))] <- fitted(obj)
  return(list(out=out,obj=obj))
}

### plot predictions
plot_arima_out <- function(out){
  plot(out[1,],ylim=range(out),type="l")
  lines(out[2,],lty=3)
  lines(out[3,],col=2)
  lines(out[4,],col=2)
  lines(out[5,],col=3)
  lines(out[6,],col=3)
}

plot_arima_obj <- function(obj){
  plot(forecast(obj,level = c(50,80),h = 25))
  lines(fitted(obj),lty=3)
  legend("topright",c("Observed","Fitted"),lty=c(1,3))
}


############################################
Method11_Saboia1977.R <- function(ASFR,
                                  joy = 1985, 
                                  obs = 30,
                                  age1 = 15, 
                                  age2 = 44,
                                  parameter = c(5,2,2),
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
  
  # extract input data: exposure
  #expos.raw <- EXPOS[paste(age1:age2),paste(year1:year2)]
  # extract input data: asfr
  fert.raw <- ASFR[paste(age1:age2),paste(year1:year2)]
  fert.raw[fert.raw<=0.00001] <- 0.00001
  raw.data <- matrix(NA,nrow=nrow(fert.raw),ncol=(year3-year1+1))
  colnames(raw.data) <- paste(year1 + c(1:ncol(raw.data))-1)
  rownames(raw.data) <- rownames(fert.raw)
  raw.data[,colnames(fert.raw)] <- fert.raw
  
  ########
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))])) ))
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
                predTFR=NA,
                predCTFR=NA,
                predCTFRM2=NA,
                predASFR=NA,
                predCASFR=NA,
                predCASFRlowerPI80 = NA,
                predASFRlowerPI80  = NA,
                predCASFRupperPI80 = NA,
                predASFRupperPI80  = NA, 
                predCASFRlowerPI95 = NA,
                predASFRlowerPI95  = NA, 
                predCASFRupperPI95 = NA,
                predASFRupperPI95  = NA,
                method="Saboia (1977)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("Saboia",paste(parameter,collapse = "_"),sep="_"),"Saboia"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  ########
  mab <- apply(fert.raw * (age1:age2 + 0.5), 2, sum) / apply(fert.raw, 2, sum)
  
  ########
  # raw tfr 
  tfr.raw <- apply(fert.raw,2,sum)
  tfr.raw <- ts(tfr.raw, start=c(year1, 1), end=c(year2, 1), frequency=1) 
  
  # Fit ARIMA
  tfr.fit <- tfr_arima_pred(tfr.raw,order=parameter,len=len,auto.bound=T)
  
  predTFR <- matrix(NA,nrow=7,ncol=(year3-year1+1))
  predTFR[1,] <- c(year1:year2,1:len+year2)
  predTFR[2:7,] <- tfr.fit$out
  colnames(predTFR) <-  paste(c(year1:year2,1:len+year2))
  rownames(predTFR) <- c("YEAR","TFR_ObsPre","TFR_FitPre","lowerPI80","upperPI80","lowerPI95","upperPI95")
  ####
  predCTFR <- predTFR
  rownames(predCTFR)[1] <- "COHORT"
  predCTFR[1,] <- predCTFR[1,] - round(mean( mab[as.character(year2 - (4:0))],na.rm = T))
  ind <- (year2-age2):(year2-age1-1)
  predCTFR[2:7,!predCTFR[1,] %in% ind] <- NA
  
  #### Set observed part as NA
  predCTFRM2 <- predCTFR
  predCTFRM2[2:7,paste(year1:year2)] <- NA
  ####
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predTFR=predTFR,
              predCTFR=predCTFR,
              predCTFRM2= predCTFRM2,
              predASFR=NA,
              predCASFR=NA,
              predCASFRlowerPI80 = NA,
              predASFRlowerPI80  = NA,
              predCASFRupperPI80 = NA,
              predASFRupperPI80  = NA, 
              predCASFRlowerPI95 = NA,
              predASFRlowerPI95  = NA, 
              predCASFRupperPI95 = NA,
              predASFRupperPI95  = NA,
              method="Saboia (1977)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("Saboia",paste(parameter,collapse = "_"),sep="_"),"Saboia"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}

