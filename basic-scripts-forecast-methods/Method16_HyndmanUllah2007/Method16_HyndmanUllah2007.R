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
#    2016.09.23
#    
#    Method: Hyndman and Ullah 2007
#     
#    Function.1: Create demogdata object from raw data matrices
#    Method16_HyndmanUllah2007.R
#    
############################################################
#### Input file: 
#     ASFR: age-specific fertility rate
#     EXPOS: age-specific exposures
#     
#### Parameters
#      joy: jump of year
#      obs: length of observation
#      age1 : star age
#      age2 : end age
#      parameter : sm.obs.var, fit.order, fit.method, fcast.method and fcast.adjust
#             1. sm.obs.var: c("theoretical") 
#                        "popconst":set the populations in every age groups are the same, and used dxt method to adjust.
#             2. fit.order: 3
#             3. fit.method: "classical"
#             4. fcast.method: "arima"
#             5. fcast.adjust: "adj"
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
#  1. Hyndman RJ, Ullah MS (2007) Robust forecasting of mortality and fertility rates: A functional
#     data approach. Computational Statistics & Data Analysis 51(10):4942-4956.
#        
#        https://robjhyndman.com/software/demography
###############################################################################
library(demography)
source("./basic-scripts-forecast-methods/Method16_HyndmanUllah2007/R.Create.demogdata.object.r")
Method16_HyndmanUllah2007.R <- function(ASFR,
                                        EXPOS,
                                        joy = 1985, 
                                        obs = 30,
                                        age1 = 15, 
                                        age2 = 44,
                                        parameter = c("theoretical","3","classical","arima","adj"),
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
  expos.raw <- EXPOS[paste(age1:age2),paste(year1:year2)]
  # extract input data: asfr
  fert.raw <- ASFR[paste(age1:age2),paste(year1:year2)]
  fert.raw[fert.raw<=0.00001] <- 0.00001
  raw.data <- matrix(NA,nrow=nrow(fert.raw),ncol=(year3-year1+1))
  colnames(raw.data) <- paste(year1 + c(1:ncol(raw.data))-1)
  rownames(raw.data) <- rownames(fert.raw)
  raw.data[,colnames(fert.raw)] <- fert.raw
  
  ########
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))]))) | !any(apply(raw.data,2,function(x) all(!is.na(x)))) )
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
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
                method="Hyndman and Ullah (2007)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("HyndmanUllah",paste(parameter,collapse = "_"),sep="_"),"HyndmanUllah"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  # generate demography object
  fert.data <- R.Create.demogdata.object.r(fert.raw,expos.raw,year1=year1,year2=year2,age1=age1,age2=age2,label=pop) 
  
  # Smooth data
  fert.sm <- smooth.demogdata(fert.data,obs.var=parameter[1])
  
  # Fit model and forecasting using ARIMA
  fert.fit <- fdm(fert.sm,order=as.numeric(parameter[2]),method=parameter[3])
  fert.fcast1 <- forecast(fert.fit,h=year3-year2,level=80,method=parameter[4],adjust=parameter[5]=="adj")
  fert.fcast2 <- forecast(fert.fit,h=year3-year2,level=95,method=parameter[4],adjust=parameter[5]=="adj")
  
 
  data <- cbind(fert.raw,fert.fcast1$rate$total)
  colnames(data) <- paste(c(as.numeric(colnames(fert.raw)),fert.fcast1$year))
  dataupper80 <- cbind(fert.raw,fert.fcast1$rate$upper)
  colnames(dataupper80) <- colnames(data)
  datalower80 <- cbind(fert.raw,fert.fcast1$rate$lower)
  colnames(datalower80) <- colnames(data)
  
  dataupper95 <- cbind(fert.raw,fert.fcast2$rate$upper)
  colnames(dataupper95) <- colnames(data)
  datalower95 <- cbind(fert.raw,fert.fcast2$rate$lower)
  colnames(datalower95) <- colnames(data)
  ####
  ################ Output step
  predASFR <- data
  predCASFR <- asfr_period_to_cohort(predASFR) 
  
  predASFRlowerPI80 <- datalower80 
  predASFRupperPI80 <- dataupper80
  predCASFRlowerPI80 <- asfr_period_to_cohort(predASFRlowerPI80)
  predCASFRupperPI80 <- asfr_period_to_cohort(predASFRupperPI80)
  
  predASFRlowerPI95 <- datalower95 
  predASFRupperPI95 <- dataupper95
  predCASFRlowerPI95 <- asfr_period_to_cohort(predASFRlowerPI95)
  predCASFRupperPI95 <- asfr_period_to_cohort(predASFRupperPI95)
  
  ################ Output step
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              predCASFRlowerPI80 = predCASFRlowerPI80,
              predASFRlowerPI80  = predASFRlowerPI80,
              predCASFRupperPI80 = predCASFRupperPI80,
              predASFRupperPI80  = predASFRupperPI80, 
              predCASFRlowerPI95 = predCASFRlowerPI95,
              predASFRlowerPI95  = predASFRlowerPI95, 
              predCASFRupperPI95 = predCASFRupperPI95,
              predASFRupperPI95  = predASFRupperPI95,
              method="Hyndman and Ullah (2007)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("HyndmanUllah",paste(parameter,collapse = "_"),sep="_"),"HyndmanUllah"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}


##############################
