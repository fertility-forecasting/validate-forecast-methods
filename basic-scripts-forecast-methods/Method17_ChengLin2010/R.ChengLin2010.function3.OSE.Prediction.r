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
#    Function.3: Estimation out-of-sample period effects and predict
#    R.ChengLin2010.function3.OSE.Prediction.r
############################################################
#### Input file: 
#    obj.APC 
#    [output of APC framework (R.ChengLin2010.function1.APC.r) or raw data]
#    obj.KPS
#    [output of KP + kernel smoothing method (R.ChengLin2010.function2.KP.smooth.r)]
#### Parameters:
#    ctfr = c("raw","bf","bfs","kp","kps")
#        raw: raw tfr
#        bf: Bongaarts and Feeney 1998
#        bfs: Bongaarts and Feeney 1998 + kernel smoothing
#        kp: Kohler and Philipov 2001
#        kps: Cheng and Lin 2010 (KP  + kernel smoothing)
#    locus.function = c("linear", "logarithmic", "quadratic")
#        patern of out-of-sample period effects
#    TX = 15
#        paramter for locus.function = "quadratic"
#### Output: obj.OSE 
#    obj.OSE$output: the same formate as obj.APC$output with updated out-of-sample predictions
#    obj.OSE$pred: prediction of ASFR (cohort * period)
#    obj.OSE$coef: estimatino of all effects
#    obj.OSE$coef.inter: estimatino of effect
#    obj.OSE$predCV   : full model Predictions (age * cohort matrix) for cohorts (year2-age1):(year2-age2+1)
#    obj.OSE$predAll  : full model predictions (age * cohort matrix) 
#    obj.OSE$predFull : full model predictions (cohort * period matrix)
#    obj.OSE$apcAll     : All predictions of APC model (age * cohort matrix) from obj.APC 
#    obj.OSE$CTFR.est   : estimation of CTFR from obj.KPS 
#    obj.OSE$coef        : all.coef
#    obj.OSE$coef.inter  : coef.inter
#    obj.OSE$coef.age    : coef.age
#    obj.OSE$coef.cohort : coef.cohort
#    obj.OSE$coef.period : coef.period
#    obj.OSE$B1  : out-of-sample period effects for each cohort, B1 for older part
#    obj.OSE$B0  : out-of-sample period effects for each cohort, B0 for younger part
#    ......
#    obj.OSE$ose.parm: parameters of APC step and current step
###############################################################################
#############################################################
R.ChengLin2010.function3.OSE.Prediction <- function(obj.APC, 
                                    obj.KPS,
                                    ctfr = "kps",
                                    locus.function = "quadratic",
                                    TX = 15, 
                                    select = "top3"
                                    ){

  

  #print("OSE.Prediction")
  ############################
  ## data
  apcout <- obj.APC$apcout
  ## parameters
  age1 <- obj.APC$age1
  age2 <- obj.APC$age2
  year1 <- obj.APC$year1
  year2 <- obj.APC$year2
  
  ########
  predCV <- matrix(NA,nrow=length(age1:age2),ncol=length((year2-age1):(year2-age2+1)))
  if(sum(is.na(obj.APC$coef))){
    return(list(pred=NA,
                predCV=predCV,
                full=NA,
                coef=NA,
                ose.parm = paste(obj.APC$apc.parm,ctfr,locus.function,TX,select,sep=".")))
  }

  ######## CTFR options
  if(ctfr == "raw"){
    CTFR.est <- obj.KPS$CTFR.Raw
  }
  if(ctfr == "bf"){
    CTFR.est <- obj.KPS$CTFR.bf
  }
  if(ctfr == "bfs"){
    CTFR.est <- obj.KPS$CTFR.bfs
  }
  if(ctfr == "kp"){
    CTFR.est <- obj.KPS$CTFR.KP
  }
  if(ctfr == "kps"){
    CTFR.est <- obj.KPS$CTFR.KPS
  }
  
  ############################
  age <- sort(unique(apcout[,"AGE"]))
  period <- sort(unique(apcout[,"YEAR"]))
  cohort <- sort(unique(apcout[,"COHORT"]))
  
  n.age <- length(age)
  n.period <- length(period)
  n.cohort <- length(cohort)  
  ############################
  ## save ASFR (real and prediction from APC) to matrix
  temp <- matrix(0,nrow=length(cohort),ncol=length(period))
  rownames(temp) <- paste(cohort)
  colnames(temp) <- paste(period)
  
  ## save corresponding index of temp
  index <- matrix(0,nrow=length(cohort),ncol=length(period))
  rownames(index) <- paste(cohort)
  colnames(index) <- paste(period)
  
  for(i in 1:nrow(apcout)){
    if(apcout[i,"YEAR"] < year1 | apcout[i,"YEAR"] > year2){
      temp[paste(apcout[i,"COHORT"]),paste(apcout[i,"YEAR"])] <- apcout[i,"PRED"]
      index[paste(apcout[i,"COHORT"]),paste(apcout[i,"YEAR"])] <- 1
    }
    if(apcout[i,"YEAR"] >= year1 & apcout[i,"YEAR"] <= year2){
      temp[paste(apcout[i,"COHORT"]),paste(apcout[i,"YEAR"])] <- apcout[i,"DEPEND"]
      index[paste(apcout[i,"COHORT"]),paste(apcout[i,"YEAR"])] <- 0
    }
  }
  
  period.use <- as.numeric(colnames(temp))
  old <- period.use[period.use > year2]
  young <- period.use[period.use < year1]
  

  ###############################
  ######## for old part
  #### Setting locus function
  if(locus.function == "linear"){
    ## locus function: linear
    locus.f <- c(1:length(old))
    TX <- length(old)
    #plot(locus.f)
    
  }
  
  if(locus.function == "logarithmic"){
    ## locus function: logarithmic
    locus.f <- log(c(1:length(old) +1))
    TX <- length(old)
    #plot(locus.f)
  }
  
  if(locus.function == "quadratic"){
    ## locus function: quadratic
    locus.f <- c(1:length(old)) * (1 - 0.5*c(1:length(old))/TX)
    locus.f[c(1:length(old)) > TX] <- locus.f[TX]
    #plot(locus.f)
  } 
  
  #### set matrix to estimate B1
  # index matrix
  out.cohort.ind <- index[,paste(year1-1)]==0 & index[,paste(year2+1)]>0 & rownames(index) %in% paste(CTFR.est[,"COHORT"])
  out.period.ind <- as.numeric(colnames(index)) > year2
  IM <- index[out.cohort.ind, out.period.ind]
  
  # locus function (sum exp{b_0 * g(x)})
  X <- IM * matrix(rep(locus.f,nrow(IM)),nrow=nrow(IM),byrow=T)
  colnames(X) <- paste("LC",colnames(X),sep="")
  
  # baseline (sum of out-of-sample prediction)
  Y <- temp[out.cohort.ind, out.period.ind]
  colnames(Y) <- paste("PD",colnames(Y),sep="")
  
  # CTFR - observed
  Z <- CTFR.est[rownames(X),"CTFR"] - apply(temp[out.cohort.ind,paste(c(year1:year2))],1,sum) + obj.APC$CONST*(age2-age1+1)

  # data set
  data <- as.data.frame(cbind(X,Y,Z))
  
  ######## use each cohort to estimate latter out-of-sample period effects (B1)
  #### nls package estimate b1: sometimes this function returns error ""Error in nlsModel() : singular""
  B1 <- NULL
  #### nlmrt package is used here
  ## nlin model
  nlsm_formula <- formula(paste0("Z ~ ", 
                                 paste(paste(paste("exp(B*",colnames(X),")",sep=""), colnames(Y), sep = "*"), collapse=" + ")
  ))
  for(i in 1:nrow(data)){
    nlsm_obj <- nlxb(nlsm_formula, data = data[i,], start = list(B=0), trace=F,lower = -5, upper = 5)
    B1 <- c(B1,summary(nlsm_obj)$coef[1])
  }
  names(B1) <- rownames(data)
  #### prediction of ASFR on old part (latter than year2)
  ind <- as.numeric(colnames(temp)) > year2
  temp[out.cohort.ind,ind] <- temp[out.cohort.ind,ind] * exp( matrix(rep(locus.f,sum(out.cohort.ind)),nrow=sum(out.cohort.ind),byrow=T) * B1) 

  ### impute cohort without CTFR estimation (older)
  out.cohort.imp.ind <- index[,paste(year1-1)]==0 & index[,paste(year2+1)]>0 & !(rownames(index) %in% paste(CTFR.est[,"COHORT"]))
  temp[out.cohort.imp.ind,ind] <- temp[out.cohort.imp.ind,ind] * exp( matrix(rep(locus.f,sum(out.cohort.imp.ind)),nrow=sum(out.cohort.imp.ind),byrow=T) * B1[length(B1)]) 
  
  ####
  old.coef <- locus.f*median(B1) + obj.APC$coef[paste("period",year2,sep="")]
  names(old.coef) <- paste("period",old,sep="")
  
  #### Impute older parts whose young parts are also imcompleted
  nlsm_obj <- nlxb(nlsm_formula, data = data[1:min(3,nrow(data)),], start = list(B=0), trace=F)
  B10 <- summary(nlsm_obj)$coef[1]
  
  out.cohort.ind <- index[,paste(year1-1)]!=0 & index[,paste(year2+1)]>0 & rownames(index) %in% paste(CTFR.est[,"COHORT"])
  if(sum(out.cohort.ind)>0){
    ind <- as.numeric(colnames(temp)) > year2
    temp[out.cohort.ind,ind] <- temp[out.cohort.ind,ind] * exp( matrix(rep(locus.f,sum(out.cohort.ind)),nrow=sum(out.cohort.ind),byrow=T) * B10) 
  }
  
  ###############################
  ######## for young part
  #### Setting locus function
  if(locus.function == "linear"){
    ## locus function: linear
    locus.f <- c(length(young):1)
    #plot(locus.f)
    
  }
  
  if(locus.function == "logarithmic"){
    ## locus function: logarithmic
    locus.f <- log(c(length(young):1 +1))
    #plot(locus.f)
  }
  
  if(locus.function == "quadratic"){
    ## locus function: quadratic
    locus.f <- c(length(young):1) * (1- 0.5*c(length(young):1)/TX)
    locus.f[c(length(young):1) > TX] <- locus.f[TX]
    #plot(locus.f)
  } 
  
  #### set matrix to estimate B0
  # index matrix
  out.cohort.ind <- index[,paste(year1-1)]>0 & rownames(index) %in% paste(CTFR.est[,"COHORT"])
  out.period.ind <- as.numeric(colnames(index)) < year1
  IM <- index[out.cohort.ind, out.period.ind]
  
  # locus function (sum exp{b_0 * g(x)})
  X <- IM * matrix(rep(locus.f,nrow(IM)),nrow=nrow(IM),byrow=T)
  colnames(X) <- paste("LC",colnames(X),sep="")
  
  # baseline (sum of out-of-sample prediction)
  Y <- temp[out.cohort.ind, out.period.ind]
  colnames(Y) <- paste("PD",colnames(Y),sep="")
  
  # CTFR - observed
  Z <- CTFR.est[rownames(X),"CTFR"] - apply(temp[out.cohort.ind,as.numeric(colnames(temp)) >= year1],1,sum)

  # data set
  data <- as.data.frame(cbind(X,Y,Z))
  
  ######## use each cohort to estimate latter out-of-sample period effects (B0)
  B0 <- NULL
  #### nls package estimate b0
  ## nlin model
  nlsm_formula <- formula(paste0("Z ~ ", 
                                 paste(paste(paste("exp(B*",colnames(X),")",sep=""), colnames(Y), sep = "*"), collapse=" + ")
  ))
  for(i in 1:nrow(data)){
    nlsm_obj <- nlxb(nlsm_formula, data = data[i,], start = list(B=0), trace=F, algorithm = "port", lower = -5, upper = 5)
    B0 <- c(B0,summary(nlsm_obj)$coef[1])
  }
  names(B0) <- rownames(data)
  
  #### prediction of ASFR on young part (earlier than year1)
  ind <- as.numeric(colnames(temp)) < year1
  temp[out.cohort.ind,ind] <- temp[out.cohort.ind,ind] * exp( matrix(rep(locus.f,nrow(data)),nrow=nrow(data),byrow=T) * B0) 

  ### impute cohort without CTFR estimation (older)
  out.cohort.imp.ind <- index[,paste(year1-1)]>0 & !(rownames(index) %in% paste(CTFR.est[,"COHORT"]))
  temp[out.cohort.imp.ind,ind] <- temp[out.cohort.imp.ind,ind] * exp( matrix(rep(locus.f,sum(out.cohort.imp.ind)),nrow=sum(out.cohort.imp.ind),byrow=T) * B0[1]) 
  
  temp <- temp - obj.APC$CONST
  temp[sign(temp) == -1] <- 0
  ###############################
  #### use "PRED" to save 
  output <- obj.APC$apcout
  ind <- output[,"YEAR"] >= year1 & output[,"YEAR"] <= year2 & output[,"AGE"] >= age1 & output[,"AGE"] <= age2
  output[!ind,"PRED"] <- NA
  for(i in 1:nrow(temp)){
    ind <- output[,"COHORT"] == as.numeric(rownames(temp))[i]
    output[ind,"PRED"] <- temp[i,paste(output[ind,"YEAR"])]
  }
  
  ####
  young.coef <- locus.f*median(B0) + obj.APC$coef[paste("period",year1,sep="")]
  names(young.coef) <- paste("period",young,sep="")
  all.coef <- c( obj.APC$coef,young.coef,old.coef)
  all.coef <- c(all.coef[1],all.coef[order(names(all.coef[-1]))+1])
  ########
  obj.APC$coef.period <- c(young.coef,obj.APC$coef.period,old.coef)
 
  ########
  ## save ASFR (real and prediction from KP2010) to matrix age*cohort
  predAll <- asfr_List_to_AC(output,vYEAR="YEAR",vAGE="AGE",vASFR="PRED")
  
  ########
  predCV <- predAll[,paste((year2-age1):(year2-age2+1))]
  
  ########
  return(list(
              predCV=predCV,
              predAll=predAll,
              predFull=temp,
              apcAll = obj.APC$apcAll,
              CTFR.est = CTFR.est,
              coef=all.coef,
              coef.inter = obj.APC$coef.inter,
              coef.age=obj.APC$coef.age,
              coef.cohort = obj.APC$coef.cohort,
              coef.period=obj.APC$coef.period,
              B1= B1,
              B0 = B0,
              age1 = obj.APC$age1, age2 = obj.APC$age2, 
              year1 = obj.APC$year1, year2 = obj.APC$year2,
              age = obj.APC$age,
              cohort = obj.APC$cohort,
              period = sort(unique(c(young,obj.APC$period,old))),
              ose.parm = paste(obj.APC$apc.parm,ctfr,locus.function,TX,select,sep=".")))
  
}

########
if(0>1){
  
  ctfr = "kps"
  locus.function = "quadratic"
  TX = 15 
  select = "top3"
  
}
