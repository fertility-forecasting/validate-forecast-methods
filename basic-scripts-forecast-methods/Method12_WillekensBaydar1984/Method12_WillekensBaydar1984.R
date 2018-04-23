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
#    2016.10.28
#    
#    Method: Willekens and Baydar 1984     
#    Method12_WillekensBaydar1984.R
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
#      parameter : NA
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
#      ASFR is adjusted by adding a small number (named as "DEPEND") 
#      In following steps, KPsmoothing and out-of-sample effects estimating, all based on "DEPEND"
#      At the last step of prediction,  this CONST will be removed
#### NOTE:
#
#### Ref:
#   1. Willekens F, Baydar N (1984) Age-period-cohort models for forecasting fertility, Working paper
#      no.45, Voorburg, July 1984.    
#        
###############################################################################
#install.packages("nlmrt")
require(nlmrt)

####
Method12_WillekensBaydar1984.R <- function(ASFR,
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
  fert.raw[fert.raw<=0.00001] <- 0.00001

  #print(c(year1,year2,age1,age2))
  #print(pop)
  ########
  if(any(is.na(sum(fert.raw[,paste(c(year1,year2))]))) | !any(apply(fert.raw,2,function(x) all(!is.na(x)))) )
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
                method="Willekens and Baydar (1984)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("WillekensBaydar",paste(parameter,collapse = "_"),sep="_"),"WillekensBaydar"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
  ################ Model step
  years <- as.numeric(colnames(fert.raw))
  ages <- as.numeric(rownames(fert.raw))
  raw.data <- expand.grid(ages,years)
  raw.data <- cbind(raw.data,as.vector(fert.raw))
  colnames(raw.data) <- c("AGE","YEAR","ASFR")

  #### sample data
  ind <- raw.data[,"YEAR"] >= year1 & raw.data[,"YEAR"] <= year2 & raw.data[,"AGE"] >= age1 & raw.data[,"AGE"] <= age2
  data <- raw.data[ind,]
  Intercept <- 1
  if(max(data[,"ASFR"]) > 50){
    CONST <- 0.001
    DEPENDlog <- log(data[,"ASFR"] + CONST)
    DEPEND <- (data[,"ASFR"] + CONST)
  }
  if(max(data[,"ASFR"]) < 20){
    CONST <- 0.000001
    DEPENDlog <- log(data[,"ASFR"] + CONST)
    DEPEND <- (data[,"ASFR"] + CONST)
  }
  COHORT <- data[,"YEAR"] - data[,"AGE"]
  data <- cbind(Intercept,data,COHORT,DEPEND,DEPENDlog)
  data <- data[order(data[,"COHORT"], data[,"AGE"]),]
  

  #### full data (used later)
  full <- expand.grid(c((year1-age2) : (year2-age1)),c(age1:age2))
  colnames(full) <- c("COHORT","AGE")
  full <- full[order(full[,"COHORT"], full[,"AGE"]),]
  
  #### Merger data and full data (impute out-sample-period)
  use <- merge(full,data,by=c("COHORT","AGE"),all.x=T,all.y=T)
  use[,"YEAR"] <- use[,"COHORT"] + use[,"AGE"]
  
  use[,"Intercept"] <- 1
  ind <- use[,"YEAR"] < year1
  use[ind, "YEAR"] <- year1
  
  ind.na <- !is.na(use[,"ASFR"])

  ######## APC function
  #### for ID-1
  ## generate dummy factors
  age <- as.factor(use[,"AGE"])
  dummy.age <- model.matrix(~0+age)
  
  
  period <- as.factor(use[,"YEAR"])
  dummy.period <- model.matrix(~0+period)
  
  for(i in 3:ncol(dummy.period)){
    dummy.period[,i] <- dummy.period[,i] + (i-2)*dummy.period[,1] - (i-1)*dummy.period[,2]
  }
  
  cohort <- as.factor(use[,"COHORT"])
  dummy.cohort <- model.matrix(~0+cohort)
  
  
  ## list of factors
  age <- sort(unique(use[,"AGE"]))
  period <- sort(unique(use[,"YEAR"]))
  cohort <- sort(unique(use[,"COHORT"]))
  
  n.age <- length(age)
  n.period <- length(period)
  n.cohort <- length(cohort)
  
  ## data for nonlinear least-square model
  data.nlin <- cbind(dummy.age,dummy.period,dummy.cohort,use)
  
  #### compute initial start value of factors based on linear model
  ## linear model
  lm_formula <- formula(paste0("DEPENDlog ~ 0+", 
                               paste(c("Intercept",
                                       colnames(dummy.age)[-1], 
                                       colnames(dummy.cohort)[-c(n.cohort)], 
                                       colnames(dummy.period)[-c(1,2)]), 
                                       collapse=" + ")
  ))
  #print("Compute initial value (ID-1)---")
  lm_obj_ID1 <- lm(lm_formula, data = data.nlin[ind.na,])
  start.value <- as.list((summary(lm_obj_ID1)$coef[,1]))
  var.name <- names(start.value)
  names(start.value) <- paste("Es_",var.name,sep="")
  
  #### Fit nlin model
  #### nls package estimate b1: sometimes this function returns error ""Error in nlsModel() : singular""
  #### nlmrt package is used here
  nlsm_formula <- formula(paste0("DEPEND ~ exp(", 
                                 paste(paste(paste("Es_",var.name,sep=""), var.name, sep = "*"), collapse=" + ")
                                 ,")"
  ))
  #print("Fit nlin model (ID-1)---")
  nlsm_obj <- nlxb(nlsm_formula, data = as.data.frame(data.nlin[ind.na,]), start = start.value, trace=F)

  coef.all.ID1 <- as.matrix(summary(nlsm_obj)$coef)
  coef.inter <- coef.all.ID1[1,1]
  coef.age <- c(0,coef.all.ID1[2:n.age,1])
  coef.cohort <- c(coef.all.ID1[n.age+1:(n.cohort-1),1],0)
  coef.period <- c(0,0,coef.all.ID1[n.age+n.cohort + 0:(year2-year1+1-3),1])
  b = -sum(coef.period * c(1:length(coef.period))) + sum(coef.period)
  a = -2*sum(coef.period) + sum(coef.period * c(1:length(coef.period)))
  coef.period <- c(a,b,coef.all.ID1[n.age+n.cohort + 0:(year2-year1+1-3),1])
  obj = auto.arima(as.ts(coef.period),max.p=5,max.d=2,max.q=2,max.P=0,max.D=0,max.Q=0)
  coef.period <- c(coef.period,forecast(obj,h=age2-age1)$mean)
  
  names(coef.inter) <- "Intercept"
  names(coef.age) <- paste("age",age,sep="")
  names(coef.cohort) <- paste("cohort",cohort,sep="")
  names(coef.period) <- paste("period",year1:(year2+age2-age1),sep="")
 
  COEF <- c(coef.inter,coef.age,coef.cohort,coef.period)
  
  #### prediction
  #PRED <- (predict(nlsm_obj,newdata = data.nlin))
  PRED <- exp(as.matrix(data.nlin[,names(COEF)]) %*% (COEF))
  RES  <- data.nlin[,"DEPEND"] - PRED
  out <- cbind(use,PRED,RES)
  ind <- !is.na(out[,"ASFR"])
  out[ind,"PRED"] <- out[ind,"ASFR"]

  ###### save ASFR (real and prediction from APC) to matrix
  ## list of factors
  age <- sort(unique(out[,"AGE"]))
  period <- sort(unique(out[,"YEAR"]))
  cohort <- sort(unique(out[,"COHORT"]))
  
  apcAll <- matrix(NA,nrow=length(age),ncol=length(cohort))
  rownames(apcAll) <- paste(age)
  colnames(apcAll) <- paste(cohort)
  
  for(i in 1:nrow(out)){
    apcAll[paste(out[i,"AGE"]),paste(out[i,"COHORT"])] <- out[i,"PRED"] - CONST
  }
  
  #### For some forecast situations, there might be extremely large values for ASFR. 
  ## The following line would set such implausible values to NA. 
  ## For the evaluation we simply use all forecasts of ASFR. 
  #apcAll[apcAll>1] <- NA
    
  ################ Output step
  predCASFR <- apcAll
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
              method="Willekens and Baydar (1984)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("WillekensBaydar",paste(parameter,collapse = "_"),sep="_"),"WillekensBaydar"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}

