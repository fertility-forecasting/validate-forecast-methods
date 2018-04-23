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
#    2017.06.07
#    
#    Method: Cheng and Lin 2010
#     
#    Function.1:  
#    Function.2:  
#    Function.3:  
#     
#    Method17_ChengLin2010.R
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
#      parameter : ctfr, locus.function
#             1. ctfr = c("raw","bf","bfs","kp","kps")
#                raw: raw tfr
#                bf: Bongaarts and Feeney 1998
#                bfs: Bongaarts and Feeney 1998 + kernel smoothing
#                kp: Kohler and Philipov 2001
#                kps: Cheng and Lin 2010 (KP  + kernel smoothing)
#             2. locus.function = c("linear", "logarithmic", "quadratic")
#                patern of out-of-sample period effects
#                paramter for locus.function = "quadratic", TX = 15
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
#
#### NOTE:
#
#### Ref:
#  1. Cheng PR, Lin ES (2010) Completing incomplete cohort fertility schedules. Demographic Research 23(9):223-256.      
#       
###############################################################################
require(nlmrt)
source("basic-scripts-forecast-methods/Method17_ChengLin2010/R.ChengLin2010.function1.APC.r")
source("basic-scripts-forecast-methods/Method17_ChengLin2010/R.ChengLin2010.function2.KP.smooth.r")
source("basic-scripts-forecast-methods/Method17_ChengLin2010/R.ChengLin2010.function3.OSE.Prediction.r")


Method17_ChengLin2010.R <- function(ASFR,
                                    joy = 1985, 
                                    obs = 30,
                                    age1 = 15, 
                                    age2 = 44,
                                    parameter = c("kp","quadratic"),
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
  # extract input data: asfr
  fert.raw <- ASFR[paste(age1:age2),paste(year1:year2)]
  #fert.raw[fert.raw<=0.00001] <- 0.00001
  raw.data <- matrix(NA,nrow=nrow(fert.raw),ncol=(year3-year1+1))
  colnames(raw.data) <- paste(year1 + c(1:ncol(raw.data))-1)
  rownames(raw.data) <- rownames(fert.raw)
  raw.data[,colnames(fert.raw)] <- fert.raw
  #print(c(year1,year2,age1,age2))
  #print(pop)
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
                method="Cheng and Lin (2010)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("ChengLin",paste(parameter,collapse = "_"),sep="_"),"ChengLin"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step

  obj.APC <- R.ChengLin2010.function1.APC(fert.raw,year1=year1,year2=year2,age1=age1,age2=age2,flag=pop,fig=F) 
  obj.KPS <- R.ChengLin2010.function2.KP.smooth(obj.APC)
  
  #KPsmooth:
  obj.pred.KPsmooth <- R.ChengLin2010.function3.OSE.Prediction(obj.APC,
                                                               obj.KPS,
                                                               ctfr = parameter[1],
                                                               locus.function = parameter[2],
                                                               TX = 15)
  
  ################ Output step
  predCASFR <- obj.pred.KPsmooth$predAll 
  predASFR <- asfr_cohort_to_period(predCASFR)
  
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              predCASFRlowerPI80 = NA,
              predASFRlowerPI80  = NA,
              predCASFRupperPI80 = NA,
              predASFRupperPI80  = NA, 
              predCASFRlowerPI95 = NA,
              predASFRlowerPI95  = NA, 
              predCASFRupperPI95 = NA,
              predASFRupperPI95  = NA,
              method="Cheng and Lin (2010)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("ChengLin",paste(parameter,collapse = "_"),sep="_"),"ChengLin"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  
}



