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
#    2016.09.07
#    
#    Method: Cheng and Lin 2010
#
#    Function.1: APC framework
#    R.ChengLin2010.function1.APC.r
############################################################
#### Input file: raw.data (with columns "AGE","YEAR","ASFR")
#### Parameters (need standard names for all methods)
#      year1: start year of sample data 
#      year2: end year of sample data
#      age1 : star age
#      age2 : end age
#      fig  : draw fig of factor effects
#      flag : name of objects and fig
#### Output: obj.APC
#      obj.APC$apcout: Raw data + prediction of ASFR (PRED) and residuals (RES) + "DEPEND"
#      obj.APC$apcCV : Predictions of APC model (age * cohort matrix) for cohorts (year2-age2):(year2-age1-1)
#      obj.APC$apcAll: All predictions of APC model (age * cohort matrix)
#      obj.APC$coef  : Intercept, age-, period- and cohort-effects
#      obj.APC$coef...  : Seperate effects
#      obj.APC$CONST    : small number used for "DEPEND"
#      obj.APC$apc.parm : parm used in this step
#### NOTE:
#      ASFR is adjusted by adding a small number (named as "DEPEND") 
#      In following steps, KPsmoothing and out-of-sample effects estimating, all based on "DEPEND"
#      At the last step of prediction,  this CONST will be removed
###############################################################################
library(nlmrt)
#####################################
asfr_AP_TabletoList <- function(raw.data){
  age <- rownames(raw.data)
  period <- colnames(raw.data)
  out <- as.data.frame(matrix(NA,nrow=length(age)*length(period),3))
  colnames(out) <- c("AGE","YEAR","ASFR")
  k=1
  for(i in 1:length(age)){
    for(j in 1:length(period)){
      out[k,] <- as.numeric(c(age[i],period[j],raw.data[i,j]))
      k <- k+1
    }
  }
  return(out)
}

asfr_List_to_AC <- function(raw.data,vYEAR="YEAR",vAGE="AGE",vASFR="ASFR"){
  COHORT <- raw.data[,vYEAR] - raw.data[,vAGE]
  temp <- cbind(raw.data[,c(vYEAR,vAGE,vASFR)],COHORT)
  colnames(temp) <- c("YEAR","AGE","ASFR","COHORT")
  
  age <- sort(unique(temp[,"AGE"]))
  period <- sort(unique(temp[,"YEAR"]))
  cohort <- sort(unique(temp[,"COHORT"]))

  ############################
  ## save ASFR (real and prediction from APC) to matrix
  out <- matrix(NA,nrow=length(age),ncol=length(cohort))
  rownames(out) <- paste(age)
  colnames(out) <- paste(cohort)
  
  for(i in 1:nrow(temp)){
    out[paste(temp[i,"AGE"]),paste(temp[i,"COHORT"])] <- temp[i,"ASFR"]
  }
  return(out)
}

asfr_AP_to_AC <- function(raw.data){
  temp <- asfr_AP_TabletoList(raw.data)
  return(asfr_List_to_AC(temp))
}
####################################
####
R.ChengLin2010.function1.APC <- function(raw.data, year1 = 1930, 
                                            year2 = 1954, 
                                            age1 = 14, 
                                            age2 = 49,
                                            fig = F,
                                            flag = "test"
                                            ){
  ########
  parm <- paste(flag,"period",year1,year2,"age",age1,age2,sep=".")
  
  ######## Data step
  if(!("ASFR" %in% colnames(raw.data))){
    raw.data <- raw.data[paste(age1:age2),paste(year1:year2)]
    if(sum(is.na(raw.data)) > 0){
      print("Non-complete data......return NA.")
      return(list(apcout=NA,
                  coef=NA,
                  coef.inter = NA,
                  coef.age=NA,
                  coef.period=NA,
                  coef.cohort = NA,
                  age=NA,period=NA,cohort=NA,
                  age1 = age1, age2 = age2, year1 = year1, year2 = year2,
                  apc.parm=parm))
    }
    age1 <- max(min(as.numeric(rownames(raw.data))),age1)
    age2 <- min(max(as.numeric(rownames(raw.data))),age2)
    raw.data <- asfr_AP_TabletoList(raw.data)
  }
  
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
  
  if(sum(is.na(data)) > 0){
    print("Non-complete data......return NA.")
    return(list(apcout=NA,
                coef=NA,
                coef.inter = NA,
                coef.age=NA,
                coef.period=NA,
                coef.cohort = NA,
                age=NA,period=NA,cohort=NA,
                age1 = age1, age2 = age2, year1 = year1, year2 = year2,
                apc.parm=parm))
  }
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
  
  ind <- use[,"YEAR"] > year2
  use[ind, "YEAR"] <- year2
  
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
  coef.period <- c(0,0,coef.all.ID1[n.age+n.cohort + 0:(n.period-3),1])
  b = -sum(coef.period * c(1:length(coef.period))) + sum(coef.period)
  a = -2*sum(coef.period) + sum(coef.period * c(1:length(coef.period)))
  coef.period <- c(a,b,coef.all.ID1[n.age+n.cohort + 0:(n.period-3),1])
  
  names(coef.inter) <- "Intercept"
  names(coef.age) <- paste("age",age,sep="")
  names(coef.cohort) <- paste("cohort",cohort,sep="")
  names(coef.period) <- paste("period",period,sep="")
  
  COEF <- c(coef.inter,coef.age,coef.cohort,coef.period)
  names(COEF) <- c("Intercept",paste("age",age,sep=""),paste("cohort",cohort,sep=""),paste("period",period,sep=""))

  #### prediction
  #PRED <- (predict(nlsm_obj,newdata = data.nlin))
  PRED <- exp(as.matrix(data.nlin[,names(COEF)]) %*% (COEF))
  RES  <- data.nlin[,"DEPEND"] - PRED
  out <- cbind(use,PRED,RES)
  out[,"YEAR"] <- out[,"COHORT"] +out[,"AGE"]
  
  ######## output figures
  ## test real values and predictions
  if(fig == T) {
    tiff(filename=paste("MethodinR/",parm,".tif",sep=""),width = 1000,height = 1000)
    
    par(mfrow=c(2,2))
    ind <- !is.na(out[,"DEPEND"])
    plot(out[ind,"DEPEND"],out[ind,"PRED"],xlab="Real ASFR",ylab="Predition of ASFR",
         main="MiR.ChengLin2010.example.APC.framework")
    abline(0,1,col="red")
    
    
    #### age effect
    plot(age,coef.age,ylim=c(-18,10),type="l",cex=2,lwd=2,col="grey",main="Age effects")
    abline(h=0)
    text(30,4,"ID-1",cex=1.5,col="grey")
    
    #### cohort effect
    plot(cohort,coef.cohort,ylim=c(-5,7),type="l",cex=2,lwd=2,col="grey",main="Cohort effects")
    abline(h=0)
    text(1910,-3,"ID-1",cex=1.5,col="grey")
    
    
    #### period effect
    plot(period,coef.period,ylim=c(-.4,.7),type="l",cex=2,lwd=2,col="grey",main="Period effects")
    abline(h=0,lwd=1.2)
    text(1940,-0.2,"ID-1",cex=1.5,col="grey")
    dev.off()
  }
  ###### save ASFR (real and prediction from APC) to matrix
  ## list of factors
  age <- sort(unique(out[,"AGE"]))
  period <- sort(unique(out[,"YEAR"]))
  cohort <- sort(unique(out[,"COHORT"]))
  
  apcAll <- matrix(NA,nrow=length(age),ncol=length(cohort))
  rownames(apcAll) <- paste(age)
  colnames(apcAll) <- paste(cohort)
  
  for(i in 1:nrow(out)){
    apcAll[paste(out[i,"AGE"]),paste(out[i,"COHORT"])] <- out[i,"PRED"]
  }
  
  apcCV <- apcAll[,paste((year2-age2):(year2-age1-1))]
    
  ####
  #print(parm)
  ####
  return(list(apcout=out,
              apcCV = apcCV,
              apcAll= apcAll,
              coef=COEF,
              coef.inter = coef.inter,
              coef.age=coef.age,
              coef.period=coef.period,
              coef.cohort = coef.cohort,
              age=age,period=period,cohort=cohort,
              age1 = age1, age2 = age2, year1 = year1, year2 = year2,
              CONST = CONST,
              apc.parm=parm))
}







#######################################################################
#### Please check this function (R.ChengLin2010.function1.APC.ID2) if needed.
#######################################################################
R.ChengLin2010.function1.APC.ID2 <- function(raw.data, year1 = 1930, 
                                         year2 = 1954, 
                                         age1 = 14, 
                                         age2 = 49,
                                         fig = F,
                                         flag = "test"
){
  
  #### compute initial start value of factors based on linear model
  ## linear model
  lm_formula <- formula(paste0("DEPENDlog ~ 0+", 
                               paste(c("Intercept",
                                       colnames(dummy.age)[-1], 
                                       colnames(dummy.period)[-c(n.period)],
                                       colnames(dummy.cohort)[-c(1,2)]), 
                                       collapse=" + ")
  ))
  print("Compute initial value (ID-2)---")
  lm_obj_ID2 <- lm(lm_formula, data = data.nlin[ind.na,])
  lm_obj_coef <- summary(lm_obj_ID2)$coef[,1]
  start.value <- as.list(lm_obj_coef)
  lower <- as.list(lm_obj_coef - 3*abs(lm_obj_coef))
  upper <- as.list(lm_obj_coef + 10*abs(lm_obj_coef))
  names(start.value) <- paste("Es_",var.name,sep="")
  names(lower) <- paste("Es_",var.name,sep="")
  names(upper) <- paste("Es_",var.name,sep="")
  
  #### Fit nlin model
  ## using nls package
  ## nlin model
  nlsm_formula <- formula(paste0("DEPEND ~ exp(", 
                                 paste(paste(paste("Es_",var.name,sep=""), var.name, sep = "*"), collapse=" + ")
                                 ,")"
  ))
  print("Fit nlin model (ID-2)---")
  nlsm_obj <- nls(nlsm_formula, 
                  data = as.data.frame(data.nlin), 
                  start = start.value, 
                  lower = lower,
                  upper = upper,
                  algorithm = "port",
                  trace=F)
  coef.all.ID1 <- summary(nlsm_obj)$coef
  coef.inter <- coef.all.ID1[1,1]
  coef.age <- c(0,coef.all.ID1[2:n.age,1])
  coef.period <- c(coef.all.ID1[n.age+ 1:(n.period-1),1],0)
  coef.cohort <- c(0,0,coef.all.ID1[n.age+n.period+0:(n.cohort-3),1])
  b = -sum(coef.cohort * c(1:length(coef.cohort))) + sum(coef.cohort)
  a = -2*sum(coef.cohort) + sum(coef.cohort * c(1:length(coef.cohort)))
  coef.cohort <- c(a,b,coef.all.ID1[n.age+n.period+0:(n.cohort-3),1])
  
  names(coef.inter) <- "Intercept"
  names(coef.age) <- paste("age",age,sep="")
  names(coef.cohort) <- paste("cohort",cohort,sep="")
  names(coef.period) <- paste("period",period,sep="")
  
  COEF <- c(coef.inter,coef.age,coef.cohort,coef.period)
  names(COEF) <- c("Intercept",paste("age",age,sep=""),paste("cohort",cohort,sep=""),paste("period",period,sep=""))
  
}

