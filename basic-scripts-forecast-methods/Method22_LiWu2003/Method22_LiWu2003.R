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
#    Method: Li and Wu 2003 
#
#    
#    Method22_LiWu2003.R
#        Stage1: SVD framework
#        Stage2: Forecasting
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
#      parameter : NA
#               
#               
#             
#      len : length of forecasting period
#      pop : population
#### Output: 
#      obsASFR: observed ASFR, age * period
#      obsCASFR: observed ASFR, age * cohort
#      predASFR: forecasted ASFR, age * period
#      predCASFR: forecasted ASFR, age * cohort
#
#
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
#   1. Li N, Wu Z (2003) Forecasting cohort incomplete fertility: A method and an application. Population Studies 57(3):303-320.   
#      
#      
###############################################################################

####################################
#install.packages("svd")
require(svd)
require(lattice)

delta <- function(x,y,z,n){
  #x:observed data
  #y:overall mean
  #z:age effect
  #n:number of observed ages
  eff <- sum(((x - y) * z)[1:n], na.rm = T) / sum(z[1:n]^2)
  return(list(eff = eff,SSres = sum((x-y-z*eff)^2),SSreg = sum((z*eff)^2)))
}

####

Rindex.thr.n <- function(R.index){
  thr.n <- 1
  if(sum(is.na(R.index))>0)
    return(NA)
  for(i in length(R.index):1){
    if(R.index[i] < 0.6){
      thr.n <- i+1
      break
    }
  }
  return(thr.n)
}

####Main function
Method22_LiWu2003.R <- function(ASFR,
                                joy = 1985, 
                                obs = 40,
                                age1 = 15, 
                                age2 = 44,
                                parameter = NA,
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
  raw.cdata <- asfr_period_to_cohort(fert.raw)
  ########
  if(any(is.na(sum(raw.data[,paste(c(year1,year2))]))) | !any(apply(raw.cdata,2,function(x) all(!is.na(x)))) )
    return(list(pop=pop,
                obsASFR=ASFR,
                obsCASFR=asfr_period_to_cohort(ASFR),
                predASFR=NA,
                predCASFR=NA,
                R.index=NA,
                R.thr.n=NA,
                method="Li Wu (2003)",
                parameter=parameter,
                label=ifelse(all(!is.na(parameter)),paste("LiWu",parameter,sep="_"),"LiWu"),
                year=c(year1,year2,year3),
                age=c(age1,age2),
                obs=obs,
                len=len,
                cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))

  ################ Fit the model
  
  ind <- apply(raw.cdata,2,function(x){sum(is.na(x))>0})
  obs.data <- raw.cdata[,!ind]

  effect.mean <- apply(obs.data,1,mean)
  obj <- propack.svd(as.matrix(obs.data - effect.mean),neig = 1)
  effect.age <- sign(obj$u[rownames(obs.data) == "22",1]) * obj$u[,1]
  effect.age <- effect.age/abs(sum(effect.age))
  
  effect.cohort <- NULL
  for(i in 1:ncol(raw.cdata)){
    effect.cohort <- c(effect.cohort,delta(raw.cdata[,i],effect.mean,effect.age,n=nrow(raw.cdata))$eff)
  }
  
  fr <- effect.age %*% t(effect.cohort) + effect.mean
  dimnames(fr) <- dimnames(raw.cdata)
  
  out <- fr
  out[!is.na(raw.cdata)] <- raw.cdata[!is.na(raw.cdata)]
  out[out < 0.00000001] <- 0.00000001
  ################

  ######## R(n) criterio
  R.index <- NULL
  for(n in 1:nrow(obs.data)){
    R.index <- c(R.index,0)
    for(i in 1:ncol(obs.data)){
      R.index[n] <- R.index[n] + delta(obs.data[,i],effect.mean,effect.age,n)$SSres
    }
  }
  R.index <-  1-R.index/sum((obs.data-effect.mean)^2)
  R.index[R.index<0 | is.na(R.index)] <- 0

  #print(R.index)
  thr.n <- 1
  for(i in length(R.index):1){
    if(R.index[i] < 0.6){
      thr.n <- i+1
      break
    }
  }
  
  
  

  ################ Output step
  predCASFR <- out
  predASFR <- asfr_cohort_to_period(predCASFR)

  return(list(pop=pop,
              obsASFR=ASFR,
              obsCASFR=asfr_period_to_cohort(ASFR),
              predASFR=predASFR,
              predCASFR=predCASFR,
              R.index=R.index,
              R.thr.n=thr.n,
              method="Li Wu (2003)",
              parameter=parameter,
              label=ifelse(all(!is.na(parameter)),paste("LiWu",parameter,sep="_"),"LiWu"),
              year=c(year1,year2,year3),
              age=c(age1,age2),
              obs=obs,
              len=len,
              cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1))))
}

################
#### test whether there is smooth (monotonic) trend in CTFR
test.monotonic.trend <- function(da,Fig=F,main=""){
  da = da[!is.na(da)]
  n.cohort = length(da)
  Cohort = names(da)
  if(is.null(Cohort)){
    Cohort <- 1:n.cohort
  }
  temp = as.data.frame(cbind(1:n.cohort,da))
  ####
  #obj <- lm(da~V1 + I(V1^2),data=temp)
  obj <- lm(da~V1 ,data=temp)
  rsq=summary(obj)$r.square
  #summary(obj)
  pred <- predict(obj,newdata = temp,se.fit = TRUE)
  pred.up95 <- qnorm(c(0.95),pred$fit,pred$se.fit)
  pred.lw95 <- qnorm(c(0.025),pred$fit,pred$se.fit)
  
  temp <- cbind(temp,pred$fit,pred$se.fit,pred.up95,pred.lw95)
  
  #### plot
  if(Fig){
    tiff(filename=paste("LiWu2003.",main,".Trend.tif",sep=""),width = 800,height = 600)
    plot(Cohort,da,type="l",ylim=range(temp[,c(2,5,6)]),ylab="Cohort total fertility",main=main)
    points(Cohort,obj$fitted.values,cex=0.3,pch=20,col=3)
    lines(Cohort,obj$fitted.values,col=3)
    lines(Cohort,pred.up95,col=2)
    lines(Cohort,pred.lw95,col=2)
    legend("topright",c("Observed","Fitted","95CI"),col=c(1,2,3),lty=1,bty="n")
    legend("top",paste("R2=",round(rsq,2)),bty="n")
    dev.off()
  }
  
  
  
  #### test "smooth" (monotonic) trend
  monotonic.trend = T
  td<- diff(temp[,"pred.up95"])
  td[abs(td)<0.001] <- 0
  if(sum(sign(td)>=0) < length(td) & sum(sign(td)<=0) < length(td)){
    monotonic.trend = F
  }
  td<- diff(temp[,"pred.lw95"])
  td[abs(td)<0.001] <- 0
  if(sum(sign(td)>=0) < length(td) & sum(sign(td)<=0) < length(td)){
    monotonic.trend = F
  }

  if(rsq<0.3){
    monotonic.trend = F
  }
  
  return(list(out=temp,obj.lm=obj,rsq=rsq,monotonic.trend=monotonic.trend))
}

########
Test.LiWu2003.assumptions <- function(obj.LiWu,Fig=F,num=10){
  if(sum(is.na(obj.LiWu$R.thr.n))>0){
    return(c(NA,NA,NA,NA))
  }
  da = apply(obj.LiWu$raw,2,sum)
  mon = test.monotonic.trend(da,Fig=Fig,main=obj.LiWu$parm)
  
  pass <- T
  if(mon$monotonic.trend==F | obj.LiWu$R.thr.n>=num){
    pass <- F
  }
  return(c(pass,mon$rsq,mon$monotonic.trend,obj.LiWu$R.thr.n))
}