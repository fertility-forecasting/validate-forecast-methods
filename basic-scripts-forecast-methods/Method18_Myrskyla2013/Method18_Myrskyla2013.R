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
#     
#    
#    Method18_Myrskyla2013.R
############################################################
#### Input file: 
#     ASFR: age-specific fertility rate
#     
#     
#### Parameters
#      joy: jump of year (first_forecast_year)
#      obs: length of observation
#      age1 : star age
#      age2 : end age
#      parameter : parameter of the method 
#            1. [years_to_extrapolate]  
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
#   1. Myrskylä, M., Goldstein, J.R. and Cheng, Y.H.A., 2013. New cohort fertility forecasts for the developed world: 
#      rises, falls, and reversals. Population and Development Review, 39(1), pp.31-56.  
#      
###############################################################################
require(fBasics)
########################################################################################################
Myrskyla2013_function <- function(data){
  n.col <- ncol(data)
  a <- data[,n.col]
  b <- (data[,n.col] - data[,1])/(n.col-1)
  kt <- apply(data-a,2,function(x){sum(x*b)/sum(b^2)})
  
  sigma <- kt[n.col]-kt[1] / (n.col-1)
  epsilon <- kt[2:n.col]-kt[2:n.col-1] - sigma
  return(list(a=a,b=b,kt=kt,sigma=sigma,epsilon=epsilon))
}
########################################################################################################
Method18_Myrskyla2013.R <- function(asfr_period_hfd,
                          joy = 2014, 
                          obs = 5,
                          age1 = 15, 
                          age2 = 40,
                          parameter = c(5),
                          len = 30,
                          folder="results/forecasts-by-method/",
                          out=F
                          ){
  
  ################ Data step
  if(is.na(joy))
    joy <- max(as.numeric(dimnames(asfr_period_hfd)[[2]]))
  
  Country <- dimnames(asfr_period_hfd)[[3]]
  ########
  year1 <- joy - obs + 1
  year2 <- joy
  year3 <- year2 + len
  
  ################ Model step
  raw.data <- array(NA,dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  raw.datalowerPI80 <- array(NA,dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  raw.datalowerPI95 <- array(NA,dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  raw.dataupperPI80 <- array(NA,dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  raw.dataupperPI95 <- array(NA,dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  var.data <- array(0, dim=c(length(age1:age2),length(year1:year3),length(Country)),dimnames=list(paste(age1:age2),paste(year1:year3),Country))
  
  ######## fitt model
  obj_all <- list()
  for(ct in Country){
    fert.raw <- asfr_period_hfd[paste(age1:age2),paste(year1:year2),ct]
    raw.data[paste(age1:age2),paste(year1:year2),ct] <- fert.raw
    raw.datalowerPI80[paste(age1:age2),paste(year1:year2),ct] <- fert.raw
    raw.datalowerPI95[paste(age1:age2),paste(year1:year2),ct] <- fert.raw
    raw.dataupperPI80[paste(age1:age2),paste(year1:year2),ct] <- fert.raw
    raw.dataupperPI95[paste(age1:age2),paste(year1:year2),ct] <- fert.raw
    obj <- Myrskyla2013_function(fert.raw[,paste((year2-4):year2)])
    obj_all[[ct]] <- obj
  }
  
  ######## estimate variance
  epsilon_Variance <- var(do.call("c",lapply(obj_all, "[[", "epsilon")),na.rm=T)
  sigma_Variance = epsilon_Variance / (5 - 1)
  
  ################ prediction step
  ######## prediction mean
  years_to_extrapolate <- as.numeric(parameter)
  if(years_to_extrapolate <=0)
    years_to_extrapolate <- len
  delta <- 1:len
  if(years_to_extrapolate < len)
    delta[years_to_extrapolate:len] <- years_to_extrapolate
  
  for(ct in Country){
    raw.data[paste(age1:age2),paste((year2+1):year3),ct] <- obj_all[[ct]]$a + as.matrix(obj_all[[ct]]$b) %*% delta
  }
  
  ######## prediction interval
  for(ct in Country){
    var.data[paste(age1:age2),paste((year2+1):year3),ct] <- as.matrix(obj_all[[ct]]$b)^2 %*%  (delta*sigma_Variance + c(1:len)*epsilon_Variance)
    raw.datalowerPI80[paste(age1:age2),paste((year2+1):year3),ct] <- raw.data[paste(age1:age2),paste((year2+1):year3),ct] - 1.281552*sqrt(var.data[paste(age1:age2),paste((year2+1):year3),ct])
    raw.datalowerPI95[paste(age1:age2),paste((year2+1):year3),ct] <- raw.data[paste(age1:age2),paste((year2+1):year3),ct] - 1.959964*sqrt(var.data[paste(age1:age2),paste((year2+1):year3),ct])
    raw.dataupperPI80[paste(age1:age2),paste((year2+1):year3),ct] <- raw.data[paste(age1:age2),paste((year2+1):year3),ct] + 1.281552*sqrt(var.data[paste(age1:age2),paste((year2+1):year3),ct])
    raw.dataupperPI95[paste(age1:age2),paste((year2+1):year3),ct] <- raw.data[paste(age1:age2),paste((year2+1):year3),ct] + 1.959964*sqrt(var.data[paste(age1:age2),paste((year2+1):year3),ct])
  }
  raw.data[raw.data<0] <- 0
  raw.datalowerPI80[raw.datalowerPI80<0] <- 0
  raw.datalowerPI95[raw.datalowerPI95<0] <- 0
  raw.dataupperPI80[raw.dataupperPI80<0] <- 0
  raw.dataupperPI95[raw.dataupperPI95<0] <- 0
  
  ######################################
  CPMobjlist <- list()
  for(ct in Country){
  ################ Output step
  predASFR <- raw.data[,,ct]
  predCASFR <- asfr_period_to_cohort(predASFR) 
  
  predASFRlowerPI80 <- raw.datalowerPI80[,,ct]
  predASFRupperPI80 <- raw.dataupperPI80[,,ct]
  predCASFRlowerPI80 <- asfr_period_to_cohort(predASFRlowerPI80)
  predCASFRupperPI80 <- asfr_period_to_cohort(predASFRupperPI80)
  
  predASFRlowerPI95 <- raw.datalowerPI95[,,ct]
  predASFRupperPI95 <- raw.dataupperPI95[,,ct]
  predCASFRlowerPI95 <- asfr_period_to_cohort(predASFRlowerPI95)
  predCASFRupperPI95 <- asfr_period_to_cohort(predASFRupperPI95)
  
  obsASFR=asfr_period_hfd[paste(age1:age2),paste((joy-obs+1):joy),ct]
  obsCASFR=asfr_period_to_cohort(obsASFR)
  CPMobjlist[[ct]] <- list(pop=ct,
                            obsASFR=obsASFR,
                            obsCASFR=obsCASFR,
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
                            method="Myrskylä et al. (2013)",
                            parameter=parameter,
                            label=ifelse(all(!is.na(parameter)),paste("Myrskyla2013",paste(parameter,collapse = "_"),sep="_"),"Myrskyla2013"),
                            year=c(year1,year2,year3),
                            age=c(age1,age2),
                            obs=obs,
                            len=len,
                            cohort=c((year1-age2),(year1-age1),(year2-age2),(year2-age1)))
  }
  
  ################ Output step
  if(!is.na(folder)){
    write_CPMobjlist_to_file(CPMobjlist,folder,header="CVlist")
  }
  if(out){
    return(CPMobjlist)
  }
}

##########################

