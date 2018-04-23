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
#    2017.02.22
#    
#    Method:     
#    Method05_Evans1986.R
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
#      parameter : parameter of the method [original/ imputed]
#             1. original: use level and pace factors calculated from ASFR of ages 15-24
#             2. imputed: impute ASFR of ages 21-24
#      len : length of forecasting period
#      pop : population
#### Output: 
#      obsASFR: observed ASFR, age * period
#      obsCASFR: observed ASFR, age * cohort
#      predASFR: forecasted ASFR, age * period
#      predCASFR: forecasted ASFR, age * cohort
#      predCASFRlowerPI[80,90] : lower boundary of PI
#      predCASFRupperPI[80/90] : upper boundary of PI
#      CF_lowerRange: lower range of cohort total fertility rate
#      CF_upperRange: upper range of cohort total fertility rate
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
#    1. Evans MDR (1986) American Fertility Patterns: A Comparison of White and Nonwhite Cohorts Born 1903-56. 
#       Population and Development Review 12(2):267-293.
#   
###############################################################################
########################################
#### Inner functions
########################################
library(forecast)

Evans_method_function_max_error <- function(x1,x2,ages=NULL){
  if(is.null(ages))
    ages <- as.numeric(colnames(x2))
  
  xx2 <- x2
  ind <- apply(x2,1,sum)
  xx2[is.na(ind),] <- NA
  #### prediction
  Pred <- xx2

  
  ####
  for(i in paste(ages)){
    obj <- lm(xx2[,i] ~ level + pace, data=as.data.frame(x1))
    lmm <- forecast(obj,newdata = as.data.frame(x1))
    Pred[,i] <- lmm$mean[ind]
  }
  ####
  CF.obsv <- apply(xx2,1,sum)
  CF.pred <- apply(Pred,1,sum)
  
  return(max(abs(CF.pred-CF.obsv)/CF.obsv,na.rm=T))
}

Evans_method_function <- function(x1,x2,ages=NULL){
  if(is.null(ages))
    ages <- as.numeric(colnames(x2))

  #### prediction
  Pred <- x2
  LMM <- array(NA,dim=c(4,ncol=ncol(x2)),dimnames=list(c("Intercept","level","pace","R2"),colnames(x2)))
  upperPI80 <- x2
  lowerPI80 <- x2
  upperPI95 <- x2
  lowerPI95 <- x2
  
  ####
  for(i in paste(ages)){
    ind <- is.na(x2[,i])
    obj <- lm(x2[,i] ~ level + pace, data=as.data.frame(x1))
    LMM[,i] <- c(obj$coef,summary(obj)$r.squared)
    
    lmm <- forecast(obj,newdata = as.data.frame(x1))
    Pred[ind,i] <- lmm$mean[ind]
    upperPI80[ind,i] <- lmm$upper[ind,1]
    lowerPI80[ind,i] <- lmm$lower[ind,1]
    upperPI95[ind,i] <- lmm$upper[ind,2]
    lowerPI95[ind,i] <- lmm$lower[ind,2]
  }

  CF.Pred=apply(Pred,1,sum)
  CF.upperPI80=apply(upperPI80,1,sum)
  CF.lowerPI80=apply(lowerPI80,1,sum)
  CF.upperPI95=apply(upperPI95,1,sum)
  CF.lowerPI95=apply(lowerPI95,1,sum)
  
  max.error <- Evans_method_function_max_error(x1,x2,ages=ages)
  wt <- apply(x2[,ages],1,function(x){sum(is.na(x))})
  max.wt <- max(wt)
  wt[which(wt==max.wt)] <- wt[which(wt==max.wt)] + c(1:sum(wt==max.wt)) - 1
  wt <- wt/max.wt
  CF.upperRange = CF.Pred * (1 + max.error * wt)
  CF.lowerRange = CF.Pred * (1 - max.error * wt)
  
  return(list(Pred=Pred,
              LMM=LMM,
              upperPI80=upperPI80,
              lowerPI80=lowerPI80,
              upperPI95=upperPI95,
              lowerPI95=lowerPI95,
              CF.Pred=CF.Pred,
              CF.upperPI80=CF.upperPI80,
              CF.lowerPI80=CF.lowerPI80,
              CF.upperPI95=CF.upperPI95,
              CF.lowerPI95=CF.lowerPI95,
              CF.upperRange =CF.upperRange,
              CF.lowerRange=CF.lowerRange
              ))
}

########################################
# Completion of cohort fertility
########################################
Method05_Evans1986.R <- function(ASFR,
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
  
  # age*period to age*cohort 
  raw.cdata <- asfr_period_to_cohort(raw.data)
  ########
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))]))) | !any(apply(raw.cdata,2,function(x) all(!is.na(x)))) )
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
              CF_lowerRange = NA,
              CF_upperRange = NA,
              method="Evans (1986)",
              parameter=parameter,
              label=ifelse(!is.na(parameter),paste("Evans",paste(parameter,collapse = "_"),sep="_"),"Evans"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  ################ Model step
  ind <- apply(raw.cdata,2,function(x){sum(is.na(x))==0})
  cohortcp <- names(ind)[ind]
  cohortinc <- min(as.numeric(cohortcp))
  obs.data <- t(raw.cdata[,paste(cohortinc : (year2-20))])  

  #### ASFRs of age >= 40 are frozen as the most recent cohort data
  if(sum(as.numeric(colnames(obs.data))>=40)>1){
    #print(obs.data)
    obs.data[,(as.numeric(colnames(obs.data))>=40)] <-  apply(obs.data[,(as.numeric(colnames(obs.data))>=40)],2,function(x){if(sum(is.na(x))>0) x[is.na(x)] <- x[max(which(!is.na(x)))]; return(x)})
  }
  
  if(parameter[1]=="imputed"){
    obs.data[,(as.numeric(colnames(obs.data))<=24)] <-  apply(obs.data[,(as.numeric(colnames(obs.data))<=24)],2,function(x){if(sum(is.na(x))>0) x[is.na(x)] <- x[max(which(!is.na(x)))]; return(x)})
  }

  #### compute level and pace for age 15-24
  cum.data <- t(apply(obs.data,1,cumsum))
  x1 <- cbind(cum.data[,"24"],cum.data[,"19"]/(cum.data[,"24"]-cum.data[,"19"]))
  colnames(x1) <- c("level","pace")

  Evans_obj <- Evans_method_function(x1,obs.data,ages=paste(25:39))

  ################ Output step
  predCASFR <- t(Evans_obj$Pred)
  predASFR <- asfr_cohort_to_period(predCASFR)
  
  predCASFRlowerPI80 <- t(Evans_obj$lowerPI80)
  predASFRlowerPI80 <- asfr_cohort_to_period(predCASFRlowerPI80) 
  predCASFRupperPI80 <- t(Evans_obj$upperPI80)
  predASFRupperPI80 <- asfr_cohort_to_period(predCASFRupperPI80)
  
  predCASFRlowerPI95 <- t(Evans_obj$lowerPI95)
  predASFRlowerPI95 <- asfr_cohort_to_period(predCASFRlowerPI95)
  predCASFRupperPI95 <- t(Evans_obj$upperPI95)
  predASFRupperPI95 <- asfr_cohort_to_period(predCASFRupperPI95)
  
  
  #### need to becarefull might get minues values:
  predCASFR[predCASFR<0] <- 0
  predASFR[predASFR<0] <- 0
  predCASFRlowerPI80[predCASFRlowerPI80<0] <- 0
  predASFRlowerPI80[predASFRlowerPI80<0] <- 0 
  predCASFRlowerPI95[predCASFRlowerPI95<0] <- 0
  predASFRlowerPI95[predASFRlowerPI95<0] <- 0
  
  predCASFRupperPI80[predCASFRupperPI80<0] <- 0
  predASFRupperPI80[predASFRupperPI80<0] <- 0
  predCASFRupperPI95[predCASFRupperPI95<0] <- 0
  predASFRupperPI95[predASFRupperPI95<0] <- 0

  #predCV <- predCASFR[,paste((year2-age1):(year2-age2+1))]
  
  ########
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
              CF_lowerRange = Evans_obj$CF.lowerRange,
              CF_upperRange = Evans_obj$CF.upperRange,
              method="Evans (1986)",
              parameter=parameter,
              label=ifelse(!is.na(parameter),paste("Evans",paste(parameter,collapse = "_"),sep="_"),"Evans"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
}

########################################################################################################


########################################################################################################






