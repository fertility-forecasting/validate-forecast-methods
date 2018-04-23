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
#    2017.02.01
#    
#    Method:     
#    Method07_Schmertmann2003.R
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
#      parameter : parameter of the method [arima/freezed]
#              arima: max ordre of (p,d,q) is (5,2,2) using auto.arima
#              freeze: use the last observed (P,H) 
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
#   1. Schmertmann CP (2003) A system of model fertility schedules with graphically intuitive pa-
#      rameters. Demographic Research 9(5):81-110.
#      http://www.demographic-research.org/Volumes/Vol12/5/
#      
###############################################################################
library(forecast)
########################################################################
####Schmertmann2003
########################################################################
Method07_Schmertmann2003.R <- function(ASFR, 
                                       joy = 1985, 
                                       obs = 30,
                                       age1 = 15, 
                                       age2 = 44,
                                       parameter = c("arima"),
                                       len = 30,
                                       pop=""
){
  #print(joy)
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
              method="Schmertmann (2003)",
              parameter=parameter,
              label=ifelse(!is.na(parameter),paste("Schmertmann2003",paste(parameter,collapse = "_"),sep="_"),"Schmertmann2003"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))

  ################ Model step
  age <- age1:age2
  cohort <- as.numeric(colnames(raw.cdata))
  ######## get knots
  ##
  AlphaAge <- rep(15,ncol(raw.cdata))
  AlphaVal <- raw.cdata[paste(15),]
  names(AlphaAge) <- names(AlphaVal) <- paste(cohort)
  ##
  tmp <- unlist(apply(raw.cdata[paste(20:31),],2,function(x){if(sum(is.na(x))>0) return(NA); max(which(x==max(x,na.rm = T)))}))
  PAge <- c(20:31)[tmp]
  PVal <- raw.cdata[cbind(as.character(PAge),as.character(cohort))]
  names(PAge) <- names(PVal) <- paste(cohort)
  ##
  HAge <- rep(NA,ncol(raw.cdata))
  HVal <- rep(NA,ncol(raw.cdata))
  names(HAge) <- names(HVal) <- paste(cohort)
  for(ch in colnames(raw.cdata)[!is.na(PAge)]){
    ind1 <- raw.cdata[-nrow(raw.cdata),ch] > PVal[ch]/2 
    ind2 <- raw.cdata[-1,ch] < PVal[ch]/2
    ind <- ind1 & ind2
    
    if(sum(ind,na.rm = T)>0){
      pos <- max(which(ind))
    }
    
    if(sum(ind,na.rm = T)==0 & sum(ind1,na.rm = T)>0){
      pos <- max(which(ind1))
    }
    
    if(c(15:44)[pos]>PAge[ch]+3){
      HAge[ch] <- c(15:44)[pos] + (PVal[ch]/2-raw.cdata[pos,ch])/(raw.cdata[pos+1,ch]-raw.cdata[pos,ch])
      HVal[ch] <- PVal[ch]/2
      #HAge[ch] <- c(15:44)[pos+1] 
      #HVal[ch] <- raw.cdata[pos+1,ch]
    }
  }
  #print(HAge)
  ######## projections of knots 
  #### arima model
  projection.of.knots.arima <- function(x){
    obj <- auto.arima(as.ts(x[(which(!is.na(x)))]),max.p=5,max.d = 2,max.q = 2)
    pre <- forecast(obj,h = length(x)-max(which(!is.na(x))))
    x[c(max(which(!is.na(x)))+1):length(x)] <- pre$mean
    return(x)
  }
  #### use linear model 
  projection.of.knots.lm <- function(x){
    ind <- !is.na(x)
    xx <- x[min(which(ind)):length(x)]
    da <- cbind(xx,as.numeric(names(xx)),(1:length(xx)))
    colnames(da) <- c("x","z","w")
    wt <- sort(c(rep((1:10)^2,4),rep(0,100)),decreasing = T)
    da[!is.na(da[,"x"]),"w"] <- sort(wt[1:sum(!is.na(da[,"x"]))])
    da[is.na(da[,"x"]),"w"] <- max(da[!is.na(da[,"x"]),"w"])
    
    obj <- lm(x~z,weights = w,data=as.data.frame(da))
    pre <- predict(obj,newdata = as.data.frame(da))
    x[max(which(ind)+1):length(x)] <- pre[names(x[max(which(ind)+1):length(x)])]
    return(x)
  }
  #### freezed rate
  projection.of.knots.freezed <- function(x){
    ind <- !is.na(x)
    if(sum(ind)>0)
    x[max(which(ind)+1):length(x)] <- x[max(which(ind))]
    return(x)
  }
  if(0>1){
    plot(PVal,type="l",ylim=c(0,0.3))
    lines(exp(projection.of.knots.arima(log(PVal))),col=2)
    lines(exp(projection.of.knots.arima(log(PVal))),col=3)
    lines(exp(projection.of.knots.arima(log(PVal))),col=4)
    lines(PVal)
  }
  
  if(parameter[1]=="arima"){
    PAge <- projection.of.knots.arima(PAge)
    PAge[PAge<20 & !is.na(PAge)] <- 20
    PAge[PAge>32 & !is.na(PAge)] <- 32
    PVal <- exp(projection.of.knots.arima(log(PVal)))
    HAge <- projection.of.knots.arima(HAge)
    HAge[HAge<=PAge+3 & !is.na(HAge)] <- PAge[HAge<=PAge+3 & !is.na(HAge)] +3
    HVal <- PVal/2
  }
  if(parameter[1]=="freeze"){
    PAge <- projection.of.knots.freezed(PAge)
    PAge[PAge<20 & !is.na(PAge)] <- 20
    PAge[PAge>32 & !is.na(PAge)] <- 32
    PVal <- projection.of.knots.freezed(PVal)
    HAge <- projection.of.knots.freezed(HAge)
    HAge[HAge<=PAge+3 & !is.na(HAge)] <- PAge[HAge<=PAge+3 & !is.na(HAge)] +3
    HVal <- PVal/2
  }

  ind <- is.na(HAge) & !is.na(PAge) 
  if(sum(ind)>0 & any(!is.na(HAge))){
    HAge[ind] <- mean(c(HAge[ind[-1]],HAge[c(F,ind[-length(ind)])]),na.rm = T)
  }
  if(sum(ind)>0 & all(is.na(HAge))){
    HAge[ind] <- 31
  }
  ################ Fit the model
  #alpha = 15 
  #P = 23
  #H = 31
  quadratic.splines.smoothing <- function(alpha,P,H,R){
    t0 <- alpha
    w <- min(0.75,0.25+0.025*(P-alpha))
    t1 <- (1-w)*alpha + w*P
    t2 <- P
    t3 <- (P+H)/2
    beta <- 50
    if(H+(H-P)/3 > 50){
      beta <- H+(H-P)/3
    }
    if(H+(H-P)*3 < 50){
      beta <- H+(H-P)*3
    }
    t4 <- (H+beta)/2
    tt <- c(t0,t1,t2,t3,t4)
    xx <- c(P,H,beta,P,beta)
    A <- matrix(c(c(xx[1]>tt)*(xx[1]-tt)^2,c(xx[2]>tt)*(xx[2]-tt)^2,c(xx[3]>tt)*(xx[3]-tt)^2,2*c(xx[1]>tt)*(xx[1]-tt),2*c(xx[3]>tt)*(xx[3]-tt)),byrow = T,nrow=5)
    theta <- solve(A)%*%c(1,0.5,0,0,0)
    out <- NULL
    for(ti in 15:44){
      out <- c(out,sum(c(ti > tt)*(ti-tt)^2 * theta) * R)
    }
    return(out) 
  }
  ########
  for(ch in  names(PVal)[!is.na(PVal)]){
    #print(ch)
    out <- quadratic.splines.smoothing(AlphaAge[ch],PAge[ch],HAge[ch],PVal[ch])
    raw.cdata[,ch] <- out[1:(age2-age1+1)]
  }

  ################ Output step
  predCASFR <- raw.cdata
  predASFR <- asfr_cohort_to_period(predCASFR)
  #predCV <- predCASFR[,paste((year2-age1):(year2-age2+1))]
  
  ########
  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              method="Schmertmann (2003)",
              parameter=parameter,
              label=ifelse(!is.na(parameter),paste("Schmertmann2003",paste(parameter,collapse = "_"),sep="_"),"Schmertmann2003"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}


########################################################################################################

########################################################################################################
