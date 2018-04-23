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
#    2017.02.10
#    
#    Method:  
#    Method04_Brass1974.R
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
#      parameter : tfr.method and par.method
#              tfr.method: projection method of TFR [input/Saboia/freeze]
#                  input("bayesTFR"): use outside TFR values 
#                  Saboia("arima"): max ordre of (p,d,q) is (5,2,2) using auto.arima
#                  freeze: use the last observed parameter (alpha, beta) 
#              par.method: projection method of parameters (alpha, beta) [arima/freeze]
#                  arima: max ordre of (p,d,q) is (5,2,2) using auto.arima
#                  freeze: use the last observed parameter (alpha, beta) 
#      len : length of forecasting period
#      pop : population
#      TFRinput: input projection of TFR
#             Row 1: TFR
#             Row 2: 80% lower boundary of PI TFR 
#             Row 3: 80% upper boundary of PI TFR 
#             Row 4: 95% lower boundary of PI TFR 
#             Row 5: 95% upper boundary of PI TFR 
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
#   1.Brass W (1974) Perspectives in Population Prediction: Illustrated by the Statistics of England
#     and Wales. Journal of the Royal Statistical Society. Series A (General) 137(4):532-583.
#    
###############################################################################
require(forecast)
########################################################################
Logit <- function(x,a=0,b=1){
  log((x-a)/(b-x))
}
Invers_Logit <- function(y,a=0,b=1){
  (exp(y)*b+a)/(exp(y)+1)
}
if(0>1){
  x <- runif(1000,1,10)
  plot(x,Logit(x,1,10))
  plot(x,Invers_Logit(Logit(x,1,10),1,10))
}
########################################################################
####Brass1978
########################################################################
Brass_trans <- function(x,sta,end){
  out <- log(-log(cumsum(x[sta:(end-1)]/sum(x[sta:end],na.rm = T))))
  out[is.infinite(out)] <- NA
  return(c(out,NA))
}
Invers_Brass_trans <- function(y,sta,end,TFR){
  xx <- exp(-exp(y[sta:(end-1)]))
  return(diff(c(0,xx*TFR,TFR)))
}
if(0>1){
  x <- c(0.00126,0.00358,0.00842,0.01265,0.02392,0.03164,0.04049,0.04801,0.05636,
         0.06445,0.07406,0.08308,0.086,0.09752,0.10206,0.10021,0.099,0.09107,0.0843,
         0.07598,0.06382,0.05258,0.04282,0.03358,0.02339,0.01673,0.01081,0.00704,0.00376,0.00214)
  y <- Brass_trans(x,1,30)
  xx <-Invers_Brass_trans(y,1,30,sum(x))
  plot(x,xx,xlab="Observed",ylab="Transformed back")
  abline(0,1)
  
  y <- read.table("R.Brass1978.data.txt",row.names = 1) 
  x <- Invers_Brass_trans(y[,1],1,40,5)
  names(x) <- 11:50
  write.table(x,"R.Brass1978.data.ASFR.txt")
  xx <- read.table("R.Brass1978.data.ASFR.txt",row.names = 1)
}

#### using linear model to estimate parameters
Brass_linear <- function (fx,fxBase,sta,end) 
{
  TFR <- sum(fx[sta:end],na.rm = T)
  yx <- Brass_trans(fx,sta,end)
  yxBase <- Brass_trans(fxBase,sta,end)
  
  obj <- lm(yx[sta:(end-1)] ~ yxBase[sta:(end-1)])
  alpha <- summary(obj)$coef[1,1]
  beta  <- summary(obj)$coef[2,1]
  
  #return(list(alpha=alpha,beta=beta))
  return(c(alpha,beta))
}
if(0>1){
  fx <- c(0.00126,0.00358,0.00842,0.01265,0.02392,0.03164,0.04049,0.04801,0.05636,
          0.06445,0.07406,0.08308,0.086,0.09752,0.10206,0.10021,0.099,0.09107,0.0843,
          0.07598,0.06382,0.05258,0.04282,0.03358,0.02339,0.01673,0.01081,0.00704,0.00376,0.00214)
  fxBase <- read.table("R.Brass1978.data.ASFR.txt",row.names = 1)[paste(15:44),1]
  out <- NULL
  for(i in c(5:30)){
    obj.Brass <- Brass_linear(fx,fxBase,1,i)
    out <- rbind(out,c(obj.Brass$alpha,obj.Brass$beta))
  }
  
  plot(c(5:30),out[,1])
  plot(c(5:30),out[,2])
  
  
  col <- rainbow(7)
  plot(fx,type="l")
  lines(fx,lw=2)
  for(i in c(5,10,15,20,25,30)){
    obj.Brass <- Brass_linear(fx,fxBase,1,i)
    fxx <- Invers_Brass_trans(obj.Brass$alpha+Brass_trans(fxBase,1,i)*obj.Brass$beta,1,i,sum(fx[1:i]))
    lines(fxx,col=col[i/5],lty=4,lw=2)
    
    fxxpred <- Invers_Brass_trans(obj.Brass$alpha+Brass_trans(fxBase,1,30)*obj.Brass$beta,1,30,sum(fx[1:i])*sum(fxBase[1:30])/sum(fxBase[1:i]))
    lines(fxxpred,col=col[i/5],lty=3,lw=2)
  }
  abline(v=c(5,10,15,20,25,30))
  legend("topright",c(paste("obs:",c(5,10,15,20,25,30),sep="")),col=col[c(5,10,15,20,25,30)/5],lty=4)
  
  
  
}


Method04_Brass1974.R <- function(ASFR,
                                 joy = 1995, 
                                 obs = 30,
                                 age1 = 15, 
                                 age2 = 40,
                                 parameter = c("arima","arima","CT"),
                                 len = 30,
                                 pop="",
                                 TFRinput=NA,
                                 the.temp.path=paste(c("U:/HFDTemp"),age2,sep="")
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
  raw.data <- matrix(NA,nrow=nrow(fert.raw),ncol=(year3-year1+1))
  colnames(raw.data) <- paste(year1 + c(1:ncol(raw.data))-1)
  rownames(raw.data) <- rownames(fert.raw)
  raw.data[,colnames(fert.raw)] <- fert.raw
  fxBase <- fert.raw[,1]
  if(parameter[3] == 'CT'){
    fxBase <- read.table("basic-scripts-forecast-methods/Method04_Brass1974/R.Brass1974.data.ASFR.txt",row.names = 1)[paste(age1:age2),1]
  }
  
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
                method="Brass (1974)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("Brass",paste(parameter,collapse = "_"),sep="_"),"Brass"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  ################ 
  ######## get TFR
  TFR <- matrix(NA,nrow=5,ncol=(year3-year1+1))
  colnames(TFR) <- paste(year1 + c(1:ncol(TFR))-1)
  TFR[1,colnames(fert.raw)] <- apply(fert.raw,2,sum,na.rm=T) 
  
  if(parameter[1] == "bayesTFR"){
    flag <- paste("obs",30,"_",joy,sep="")
    obj.years <- c((joy-obs+1):joy)
    out.file <- file.path(the.temp.path, flag,paste("bayesTFR.pred",flag,"Rdata",sep="."))
    out.file <- file.path(the.temp.path, flag,"predictions","prediction.rda")
    print(out.file)
    rm(bayesTFR.prediction)
    load(out.file)
    daobs <- bayesTFR.prediction$tfr_matrix_reconstructed
    dapre <- bayesTFR.prediction$quantiles[,c("0.5","0.1","0.9","0.025","0.975"),]

    
    info <- read.table(file.path(the.temp.path, "include.txt"),header=T)

    
    ct <- which(info[,"name"] == full_names[pop])

    if(!info[ct,"country_code"] %in% colnames(daobs)){
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
                  method="Brass (1974)",
                  parameter=parameter,
                  label=ifelse(all(!is.na(parameter)),paste("Brass",paste(parameter,collapse = "_"),sep="_"),"Brass"),
                  year=c(year1,year2,year3),
                  age=c(age1,age2),
                  obs=obs,
                  len=len,
                  cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
    }
    
    i <- which (colnames(daobs) == info[ct,"country_code"])
    yr <- intersect(colnames(TFR)[is.na(TFR[1,])],colnames(dapre[i,,]))
    TFR[,yr] <- dapre[i,,yr]
  }
  
  ######## get paramerters
  param <- matrix(NA,nrow=2,ncol=(year3-year1+1))
  colnames(param) <- paste(year1 + c(1:ncol(param))-1)
  param[1:2,colnames(fert.raw)] <- apply(fert.raw,2,Brass_linear,fxBase,1,age2-age1+1)
  
  ######## projections of knots 
  #### arima model including PI
  projection.of.knots.arima.PI <- function(x){
    obj <- auto.arima(as.ts(x[!is.na(x)]),max.p=5,max.d = 2,max.q = 2)
    pre <- forecast(obj,h = length(x)-max(which(!is.na(x))))
    out <- t(cbind(pre$mean,pre$lower[,1],pre$upper[,1],pre$lower[,2],pre$upper[,2]))
    #out[out < 0.00001] <- 0.00001 ## this might be used as an alternative lower boundary of the parameters
    out[out < 0.01] <- 0.01
    return(out)
  }
  #### arima model
  projection.of.knots.arima <- function(x){
    obj <- auto.arima(as.ts(x[!is.na(x)]),max.p=5,max.d = 2,max.q = 2)
    pre <- forecast(obj,h = length(x)-max(which(!is.na(x))))
    x[c(max(which(!is.na(x)))+1):length(x)] <- pre$mean
    ##x[x < 0.00001] <- 0.00001 ## this might be used as an alternative lower boundary of the parameters
    x[x < 0.01] <- 0.01 ## this was also used in the first version of CV
    return(x)
  }
  #### freezed rate
  projection.of.knots.freezed <- function(x){
    ind <- !is.na(x)
    x[max(which(ind)+1):length(x)] <- x[max(which(ind))]
    return(x)
  }
  
  ######## prejection of TFR 
  tfr.method <- parameter[1]
  if(!tfr.method %in% c("Saboia","freeze") & parameter[1] != "bayesTFR"){
    yr <- intersect(colnames(TFR),colnames(TFRinput))
    TFR[,yr] <- TFRinput[,yr]
  }
  if(tfr.method=="Saboia"){
    TFR[,c(max(which(!is.na(TFR[1,])))+1):length(TFR[1,])] <- projection.of.knots.arima.PI(TFR[1,])
  }
  if(tfr.method=="freeze"){
    TFR[1,] <- projection.of.knots.freezed( TFR[1,])
  }
  
  ######## prejection of parameter
  par.method <- parameter[2]
  if(par.method=="arima"){
    param[1,] <- projection.of.knots.arima(param[1,])
    param[2,] <- projection.of.knots.arima(param[2,])
  }
  if(par.method=="freeze"){
    param[1,] <- projection.of.knots.freezed( param[1,])
    param[2,] <- projection.of.knots.freezed( param[2,])
  }
  
  ################
  predASFR <- raw.data
  predASFRlowerPI80 <- raw.data 
  predASFRupperPI80 <- raw.data
  predASFRlowerPI95 <- raw.data 
  predASFRupperPI95 <- raw.data
  
  ##### project of profile
  for(yr in colnames(raw.data[,is.na(raw.data[1,])])){
    predASFR[,yr]          <- Invers_Brass_trans(param[1,yr]+Brass_trans(fxBase,1,age2-age1+1)*param[2,yr],1,age2-age1+1,TFR[1,yr])
    predASFRlowerPI80[,yr] <- Invers_Brass_trans(param[1,yr]+Brass_trans(fxBase,1,age2-age1+1)*param[2,yr],1,age2-age1+1,TFR[2,yr])
    predASFRupperPI80[,yr] <- Invers_Brass_trans(param[1,yr]+Brass_trans(fxBase,1,age2-age1+1)*param[2,yr],1,age2-age1+1,TFR[3,yr])
    predASFRlowerPI95[,yr] <- Invers_Brass_trans(param[1,yr]+Brass_trans(fxBase,1,age2-age1+1)*param[2,yr],1,age2-age1+1,TFR[4,yr])
    predASFRupperPI95[,yr] <- Invers_Brass_trans(param[1,yr]+Brass_trans(fxBase,1,age2-age1+1)*param[2,yr],1,age2-age1+1,TFR[5,yr])
  }
  
  ########
  predCASFR <- asfr_period_to_cohort(predASFR) 
  predCASFRlowerPI80 <- asfr_period_to_cohort(predASFRlowerPI80)
  predCASFRupperPI80 <- asfr_period_to_cohort(predASFRupperPI80)
  predCASFRlowerPI95 <- asfr_period_to_cohort(predASFRlowerPI95)
  predCASFRupperPI95 <- asfr_period_to_cohort(predASFRupperPI95)
  
  ########
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
              method="Brass (1974)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("Brass",paste(parameter,collapse = "_"),sep="_"),"Brass"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
  
}

