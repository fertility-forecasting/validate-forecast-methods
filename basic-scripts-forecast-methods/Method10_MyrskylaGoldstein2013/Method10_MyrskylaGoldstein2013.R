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
#    Method10_MyrskylaGoldstein2013.R
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
#            1. [origial/ freezed / arima] method
#            2. [T/F] IFC.adj: infecundity correction or not (T/F)
#            3. [param/asfr] prior.opt : used freezed/arima method on parameter or ASFR
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
#   1. Myrskyla, Mikko, and Joshua R. Goldstein. "Probabilistic forecasting using stochastic diffusion models, 
#      with applications to cohort processes of marriage and fertility." Demography 50.1 (2013): 237-260.
#   
###############################################################################
require(forecast)
####################################################
#### Inner functions
####################################################
## generate g_t
gt_function <- function(x,model="Hernes"){
  if(model == "Hernes")
    gt <- log(1/(2 * x[-c(1,length(x))] * (1-x[-c(1,length(x))])) * diff(x,lag=2))

  if(model == "Gompertz")
    gt <- log(1/(2 * x[-c(1,length(x))]) * diff(x,lag=2))
  
  if(model == "logistic")
    gt <- log(1/(2 * x[-c(1,length(x))]^2) * diff(x,lag=2))
  
  return(gt)
}
#gt <- gt_function(x)
####################################################
## estimation of the parameters
param_est <- function(gt){
  gt <- gt[!is.na(gt)]
  delta <- (gt[length(gt)] - gt[1])/(length(gt)-1)
  sigma2 <- sum((diff(gt)-delta)^2)/(length(gt)-3)
  return(c(delta,sigma2))
}
#param_est(gt)
####################################################
## prediction of g_t
gt_pred <- function(gt_ini,delta,len=30,ifc=1){
  gtk <- gt_ini + delta * cumsum(ifc^c(1:len))
  return(gtk)
}

####################################################
## prediction
MG_prediction <- function(x_ini,gt_ini,gt,delta,sigma2,model="Hernes"){
  len <- length(gt)
  if(model == "Hernes"){
    ## prediction of x
    xx <- c(x_ini,rep(NA,len))
    for(i in 1:len+1){
      xx[i] <- xx[i-1] + xx[i-1] * (1-xx[i-1]) * exp(gt[i-1])
    }
    ## prediction of variance
    gam <- xx * (1-xx) # gamma
  }

  if(model == "Gompertz"){
    ## prediction of x
    xx <- c(x_ini,rep(NA,len))
    for(i in 1:len+1){
      xx[i] <- xx[i-1] /(1- exp(gt[i-1])) 
    }
    ## prediction of variance
    gam <- rep(1,length(xx))
  }
  
  if(model == "logistic"){
    ## prediction of x
    xx <- c(x_ini,rep(NA,len))
    for(i in 1:len+1){
      xx[i] <- xx[i-1] + xx[i-1]^2  * exp(gt[i-1])
    }
    ## prediction of variance
    gam <- xx^2
  }
  
  #### 
  covm <- matrix(0,nrow=length(xx),ncol=length(xx))  # covariance matrix
  for(i in 1:length(xx))
    for(j in 1:length(xx)){
      covm[i,j] <- min(i,j) * exp(delta*(i+j)) * gam[i]*gam[j]
    }
  
  vx <- c(0,rep(0,len))
  for(i in 1:len+1){
    vx[i] <- sigma2 * exp(2*gt_ini) * sum(covm[1:(i-1),1:(i-1)])
  }
  ## prediction with prediction interval
  results <- matrix(NA,nrow=6,ncol=len)
  rownames(results) <- c("Pred","Var","upperPI80","lowerPI80","upperPI95","lowerPI95")
  
  results[1,] <- xx[-1]
  results[2,] <- vx[-1]
  results[3,] <- qnorm(0.80, mean = xx,sd = sqrt(vx))[-1]
  results[4,] <- qnorm(0.10, mean = xx,sd = sqrt(vx))[-1]
  results[5,] <- qnorm(0.975,mean = xx,sd = sqrt(vx))[-1]
  results[6,] <- qnorm(0.025,mean = xx,sd = sqrt(vx))[-1]
  
  return(results)

}

MG_method_one_cohort <- function(x,model="Hernes",len=30,ifc=1){
  gt <- gt_function(x,model=model)
  param <- param_est(gt)
  gt_pred <- gt_pred(gt[length(gt)],param[1],len=len,ifc=ifc)
  out <- MG_prediction(x[length(x)],gt[length(gt)],gt_pred,param[1],param[2],model=model)
  return(list(results=out,x=x,gt=gt,param=param,len=len))
}
#obj <- MG_method_one_cohort(x)


#### estimate IFC
MG_method_ifc_estimation <- function(x,model="Gompertz"){
  gt <- apply(x,1,gt_function,model)
  ages <- as.numeric(rownames(gt))
  agemax <- max(ages)
  
  param <- t(apply(as.matrix(gt[paste(ages[ages<=30]),]),2,param_est))
  
  wt <- (apply(x,1,function(x){(agemax-c(31:agemax)+1)/((agemax-30)*8)}))
  
  gt_pred_mat <- function(x,len=30,ifc=1){
    gt_pred(x[1],x[2],len=len,ifc=ifc)
  }
  
  obs_gt31to45 <- as.matrix(gt[paste(31:agemax),])
  wse <- NULL
  for(ifc in seq(1,1.5,by=0.0001)){
    pred_gt31to45 <- (apply(as.matrix(cbind(gt["30",],param)),1,gt_pred_mat,len=nrow(obs_gt31to45),ifc))
    wse <- c(wse, sum((pred_gt31to45 - obs_gt31to45)^2 *wt,na.rm=T))
  }
  ifc <- seq(1,1.5,by=0.0001)[which(wse==min(wse))]
  
  return(ifc)
}

#### prediction
#require(forecast)
MG_method_prediction <- function(x,model="Gompertz",ifc=1,method="origial",prior.opt = "param",threshold=3){
  xx <- x
  ## impute un-observed values untill age 31
  if(method == "freezed" & prior.opt %in% c("asfr","ASFR")){
    for(i in colnames(xx)){
      if(as.numeric(i) <= 31){
        ind <- is.na(xx[,i])
        xx[ind,i] <- xx[max(which(!ind)),i]
      }
    }
  }
  
  if(method == "arima" & prior.opt %in% c("asfr","ASFR")){
    for(i in colnames(xx)){
      if(as.numeric(i) <= 31){
        ind <- is.na(xx[,i])
        if(sum(ind)>0)
          xx[ind,i] <- forecast(auto.arima(as.ts(xx[!ind,i]),max.p = 5,max.d = 2,max.q = 2),h = sum(ind))$mean
      }
    }
  }
  
  gt <- cbind(t(apply(xx,1,gt_function,model)),NA)
  colnames(gt)[ncol(gt)] <- colnames(xx)[ncol(xx)]
  ages <- as.numeric(colnames(gt))
  agemax <- max(ages)
  
  param <- t(apply(gt[,paste(ages[ages<=30])],1,param_est))
  
  ## projecting parameters (delta,sigma) using auto.arima model based on oberved data
  if(method == "arima" & prior.opt == "param"){
    ind <- !is.na(gt[,"30"])
    if(model=="Gompertz"){
      obj <- forecast(auto.arima(as.ts(1/(1-exp(param[ind,1]))),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean ## CFT
      if(any(obj<0 | is.na(obj) | obj > 20)){
        param[!ind,1] <- forecast(auto.arima(as.ts(param[ind,1]),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean
      }
      if(!any(obj<0 | is.na(obj) | obj > 20)){
        param[!ind,1] <- log(max(1-1/(forecast(auto.arima(as.ts(1/(1-exp(param[ind,1]))),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean + 1e-20), 0.000001))
      }
    }
    if(model=="Hernes"){  
      param[!ind,1] <- log(max(forecast(auto.arima(as.ts(exp(param[ind,1])),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean,0.00000001))
    }
    if(model=="logistic"){ 
      param[!ind,1] <- log(forecast(auto.arima(as.ts(exp(param[ind,1])),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean)
    }
    
    param[!ind,2] <- max(exp(forecast(auto.arima(as.ts(log(param[ind,2]+0.00001)),max.p = 5,max.d = 2,max.q = 2),h = sum(!ind))$mean)-0.00001,0.000001)
  }
  
  if(method == "freezed" & prior.opt == "param"){
    ind <- !is.na(gt[,"30"])
    param[!ind,1] <- param[max(which(ind)),1]
    param[!ind,2] <- param[max(which(ind)),2]
  }
  
  #### prediction of gt
  for(r in rownames(gt)){
    agelst <- max(ages[!is.na(gt[r,])])
    if(agelst<30){
      gt[r,paste((agelst+1):30)] <- gt_pred(gt[r,paste(agelst)],param[r,1],len=30-agelst,ifc=1)
    }
    
    agelst <- max(ages[!is.na(gt[r,])])
    if(agelst < max(ages)){
      gt[r,paste((agelst+1):agemax)] <- gt_pred(gt[r,paste(agelst)],param[r,1],len=agemax-agelst,ifc=ifc)
    }
  }
  
  #### prediction of Pt
  Pred <- x
  Var <- x
  Var[!is.na(Var)] <- 0
  upperPI80 <- x
  lowerPI80 <- x
  upperPI95 <- x
  lowerPI95 <- x
  
  for(r in rownames(Pred)[is.na(Pred[,paste(agemax)])]){
    agelst <- !is.na(Pred[r,])
    agelst <- as.numeric(max(names(agelst[agelst])))

    obj <- MG_prediction(Pred[r,paste(agelst)],gt[r,paste(agelst)],gt[r,paste((agelst+1):agemax)],param[r,1],param[r,2],model=model)
    Pred[r,paste((agelst+1):agemax)] <- obj[1,]
    Var[r,paste((agelst+1):agemax)]  <- obj[2,]
    upperPI80[r,paste((agelst+1):agemax)] <- obj[3,]
    lowerPI80[r,paste((agelst+1):agemax)] <- obj[4,]
    upperPI95[r,paste((agelst+1):agemax)] <- obj[5,]
    lowerPI95[r,paste((agelst+1):agemax)] <- obj[6,]
  }
  
  coln <- colnames(Pred)
  cfr <- Pred[,ncol(Pred)]
  max_cfr <- max(cfr[!is.na(x[,"44"])])
  Pred <- cbind(Pred[,1],t(apply(Pred,1,diff)))
  upperPI80 <- cbind(upperPI80[,1],t(apply(upperPI80,1,diff)))
  lowerPI80 <- cbind(lowerPI80[,1],t(apply(lowerPI80,1,diff)))
  upperPI95 <- cbind(upperPI95[,1],t(apply(upperPI95,1,diff)))
  lowerPI95 <- cbind(lowerPI95[,1],t(apply(lowerPI95,1,diff)))
  
  ##### Set boundaries to remove outliers

  for(ag in paste(20:45)){
    max_pred <- max(Pred[!is.na(x[,ag]),ag],na.rm = T)
    
    outlier_max_ind <- (Pred[,ag]) > max_pred * threshold & cfr > min(max_cfr * threshold, 15)
    if(sum(outlier_max_ind) > 0){
      Pred[outlier_max_ind,ag] <- max_pred * threshold
      upperPI80[outlier_max_ind,ag] <- max_pred * threshold
      lowerPI80[outlier_max_ind,ag] <- max_pred * threshold
      upperPI95[outlier_max_ind,ag] <- max_pred * threshold
      lowerPI95[outlier_max_ind,ag] <- max_pred * threshold
    }
  }


  
  Pred[Pred<0] <- 0
  lowerPI80[lowerPI80<0] <- 0
  lowerPI95[lowerPI95<0] <- 0
  
  
  colnames(Pred) <- coln
  colnames(upperPI80) <- coln
  colnames(lowerPI80) <- coln
  colnames(upperPI95) <- coln
  colnames(lowerPI95) <- coln
  
  return(list(Pred=Pred,Var=Var,upperPI80=upperPI80,lowerPI80=lowerPI80,upperPI95=upperPI95,lowerPI95=lowerPI95))
}




########################################
# Completion of cohort fertility
########################################
Method10_MyrskylaGoldstein2013.R <- function(ASFR,
                                             joy = 1985, 
                                             obs = 30,
                                             age1 = 15, 
                                             age2 = 44,
                                             parameter = c("original",T,"param"),
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
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))]))))
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
                method="Myrskylä and Goldstein (2013)",
                parameter=parameter,
                label=ifelse(!is.na(parameter),paste("MG2013",paste(parameter,collapse = "_"),sep="_"),"MG2013"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1)))) 
  ################ Model step
  cum.data <- apply(raw.cdata,2,cumsum)
  ind <- apply(cum.data,2,function(x){sum(is.na(x))==0})
  if(sum(ind)==0)
    ind <- apply(cum.data[paste(max(15,age1):min(43,age2)),],2,function(x){sum(is.na(x))==0})
  cohortcp <- names(ind)[ind]

  #### infecundity correction or not
  ifc <- 1
  if(as.logical(parameter[2])){
    if(length(cohortcp)>1){
      obs.data <- t(cum.data[,cohortcp]) ### select the last 10 cohorts
      ifc <- MG_method_ifc_estimation(obs.data)
    }
    if(length(cohortcp)==1){
      obs.data <- as.matrix(t(cum.data[,cohortcp]))
      rownames(obs.data) <- paste(cohortcp)
      ifc <- MG_method_ifc_estimation(obs.data)
    }
    
  }
    
  
  cohortinc <- min(as.numeric(cohortcp)) + 1
  x <- t(cum.data[,paste(cohortinc : (year2-20))])
  MG_obj <- MG_method_prediction(x,model = "Gompertz", ifc = ifc, method=parameter[1], prior.opt=parameter[3])
  
  
  ################ Output step
  predCASFR <- t(MG_obj$Pred)
  predASFR <- asfr_cohort_to_period(t(MG_obj$Pred))

  predASFRlowerPI80 <- asfr_cohort_to_period(t(MG_obj$lowerPI80))
  predASFRupperPI80 <- asfr_cohort_to_period(t(MG_obj$upperPI80))
  predCASFRlowerPI80 <- t(MG_obj$lowerPI80)
  predCASFRupperPI80 <- t(MG_obj$upperPI80)
  
  predASFRlowerPI95 <- asfr_cohort_to_period(t(MG_obj$lowerPI95))
  predASFRupperPI95 <- asfr_cohort_to_period(t(MG_obj$upperPI95))
  predCASFRlowerPI95 <- t(MG_obj$lowerPI95)
  predCASFRupperPI95 <- t(MG_obj$upperPI95)
  
  predCASFR[predCASFR<0] <- 0
  predASFR[predASFR<0] <- 0
  predASFRlowerPI80[predASFRlowerPI80<0] <- 0
  predASFRlowerPI95[predASFRlowerPI95<0] <- 0
  
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
              method="Myrskylä and Goldstein (2013)",
              parameter=parameter,
              label=ifelse(!is.na(parameter),paste("MG2013",paste(parameter,collapse = "_"),sep="_"),"MG2013"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
}

########################################################################################################

########################################################################################################
