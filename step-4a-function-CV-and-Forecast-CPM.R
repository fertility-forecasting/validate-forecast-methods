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


#########################################################################################################################
#    CV
#
#
#
#
#
#
#
#
#
#########################################################################################################################
########
write_CPMobj_to_file <- function(CPMobj,folder=NA,pop,header="obj"){
  if(!dir.exists(folder))
    dir.create(folder)
  
  flag <- paste(header,CPMobj$label,"obs",CPMobj$obs,"joy",CPMobj$year[2],"len",CPMobj$len,CPMobj$pop,sep="_")
  save(CPMobj,file = ifelse(!is.na(folder),file.path(folder,paste(flag,".Rdata",sep="")),paste(flag,".Rdata",sep="")))
}
########
write_CPMobjlist_to_file <- function(CPMobjlist,folder=NA,header="list"){
  if(!dir.exists(folder))
    dir.create(folder)
  
  flag <- paste(header,CPMobjlist[[1]]$label,"obs",CPMobjlist[[1]]$obs,"joy",CPMobjlist[[1]]$year[2],"len",CPMobjlist[[1]]$len,sep="_")
  save(CPMobjlist,file = ifelse(!is.na(folder),file.path(folder,paste(flag,".Rdata",sep="")),paste(flag,".Rdata",sep="")))
}
########
##########################################################
asfr_period_to_cohort <- function(input){
  temp <- expand.grid(as.numeric(dimnames(input)[[1]]),as.numeric(dimnames(input)[[2]]))
  cohort <- temp[,2] - temp[,1]
  output <- tapply(as.vector(as.matrix(input)), list(temp[,1], cohort),sum)
  return(as.data.frame(output))
}

asfr_cohort_to_period <- function(input,clean=T){
  temp <- expand.grid(as.numeric(dimnames(input)[[1]]),as.numeric(dimnames(input)[[2]]))
  period <- temp[,2] + temp[,1]
  output <- tapply(as.vector(as.matrix(input)), list(temp[,1], period),sum)
  if(clean){
    output <- output[,!apply(output,2,function(x) all(is.na(x)))]
  }
  return(as.data.frame(output))
}


###########################################################
extract_obs_pred_data <- function(CVobse,CVpred,age1=15,age2=40,AAF1=20,AAF2=age2-1,flag="test",out=""){
  pop_names <- dimnames(CVobse)[[3]]
  names(pop_names) <- pop_names
  cfr_observed  <- array(NA,dim=c(length(AAF1:AAF2),length(pop_names),length(1895:2014)),dimnames=list(AAF1:AAF2,pop_names,2014:1895)) 
  cfr_forecast  <- array(NA,dim=c(length(AAF1:AAF2),length(pop_names),length(1895:2014)),dimnames=list(AAF1:AAF2,pop_names,2014:1895)) 
  
  
  for(pop in pop_names){
    #print(pop)
    temp <- CVobse[,,pop]
    period <- as.numeric(names(temp[1,])[!is.na(temp[1,])])
    ctemp <- asfr_period_to_cohort(temp)
    pred <- CVpred[pop,,,]
    years <- as.numeric(dimnames(pred)[[1]])
    for(joy in years){
      for(AAF in AAF1:AAF2){
        ####
        #print(c(AAF,joy,joy-AAF))
        cfr_observed[as.character(AAF), pop_names[pop],as.character(joy)] <- sum(ctemp[as.character(age1:age2),as.character(joy-AAF)]) 
        cfr_forecast[as.character(AAF), pop_names[pop],as.character(joy)] <- sum(pred[as.character(joy),as.character(age1:age2),as.character(joy-AAF)]) 
      } 
    }
  }
  
  
  file <- file.path(out,paste(flag,"_cfr_CVobse.Rdata",sep=""))
  save(cfr_observed,file = file)
  file <- file.path(out,paste(flag,"_cfr_CVpred.Rdata",sep=""))
  save(cfr_forecast,file = file)
}

###########################################################
save_observed_predicted_cfr <- function(asfr_period_data,obs=obs,age1=15,age2=40,parameter= c("linear","arima"),len=50,name="CoaleTrussell",file.head="Forecastlist",folder=""){
  
  label <- ifelse(all(!is.na(parameter)),paste(name,paste(parameter,collapse = "_"),sep="_"),name)
  if(name == "")
    label <- ifelse(all(!is.na(parameter)),paste(paste(parameter,collapse = "_"),sep="_"),name)
  
  flag <- paste("Forecastlist",label,"obs",obs,"len",len,sep="_")
  file <- file.path(folder,paste(flag,"_CVpredCASFR.Rdata",sep=""))
  
  if(exists("CVpredCASFR")){
    rm("CVpredCASFR")
  }
  
  if(exists("cfr_forecast")){
    rm("cfr_forecast")
  }
  load(file)
  cat("Calculate CFR...")
  ####
  extract_obs_pred_data(CVobse =asfr_period_data, CVpred = CVpredCASFR,flag="Result",age1=age1,age2=40,out=folder)
  cat("done.\n")
}

########################
# Cross-validation
########################
CV_CPM_ASFR_function <- function(FUN, asfr_period_hfd,joy=2010,obs=30,age1=15,age2=44,parameter=NA,len=50,folder=NA,out=F){

  CPMobjlist <- apply(asfr_period_hfd,3,FUN,joy,obs,age1,age2,parameter,len)
  for(pop in dimnames(asfr_period_hfd)[[3]]){
    CPMobjlist[[pop]]$pop <- pop
  }
  if(!is.na(folder)){
    write_CPMobjlist_to_file(CPMobjlist,folder,header="Forecastlist")
  }
  if(out){
    return(CPMobjlist)
  }
}

if(0>1){
  source("Scripts/Method00_FreezeRate/Method00_FreezeRate.R")
  pop = "AUT"
  ASFR = asfr_period_hfd[,,pop]
  joy = 1995 
  obs = 30
  age1 = 15 
  age2 = 44
  parameter = NA
  len = 50
  
  CV_CPM_ASFR_function(Method00_FreezeRate.R,asfr_period_hfd,joy=2010,obs=30,age1=15,age2=44,parameter=NA,len=50,folder="Temp")
}

########################
# Real-forecast
########################
Forecast_CPM_ASFR_function <- function(FUN, asfr_period_hfd,joy=NA,obs=30,age1=15,age2=44,parameter=NA,len=50,folder=NA,out=F){
  
  #CPMobjlist <- apply(asfr_period_hfd,3,FUN,joy=NA,obs,age1,age2,parameter,len)
  #for(pop in dimnames(asfr_period_hfd)[[3]]){
  #  CPMobjlist[[pop]]$pop <- pop
  #}
  CPMobjlist <- list()
  for(pop in dimnames(asfr_period_hfd)[[3]]){
    #print(pop)
    CPMobjlist[[pop]] <- FUN(asfr_period_hfd[,,pop],joy,obs,age1,age2,parameter,len,pop)
    CPMobjlist[[pop]]$pop <- pop
  }
  
  if(!is.na(folder)){
    write_CPMobjlist_to_file(CPMobjlist,folder,header="Forecastlist")
  }
  if(out){
    return(CPMobjlist)
  }
}

if(0>1){
  source("Scripts/Method00_FreezeRate/Method00_FreezeRate.R")
  pop = "AUT"
  ASFR = asfr_period_hfd[,,pop]
  joy = 1995 
  obs = 30
  age1 = 15 
  age2 = 44
  parameter = NA
  len = 50
  
  test <- Forecast_CPM_ASFR_function(Method00_FreezeRate.R,asfr_period_hfd,joy=NA,obs=30,age1=15,age2=44,parameter=NA,len=50,folder="Temp",out=T)
}

Forecast_CPM_ASFR_EXPOS_function <- function(FUN, asfr_period_hfd,expos_period_hfd,joy=NA,obs=30,age1=15,age2=44,parameter=NA,len=50,folder=NA,out=F){
  
  CPMobjlist <- list()
  for(pop in dimnames(asfr_period_hfd)[[3]]){
    CPMobjlist[[pop]] <- FUN(asfr_period_hfd[,,pop],expos_period_hfd[,,pop],joy,obs,age1,age2,parameter,len,pop)
    CPMobjlist[[pop]]$pop <- pop
  }
  if(!is.na(folder)){
    write_CPMobjlist_to_file(CPMobjlist,folder,header="Forecastlist")
  }
  if(out){
    return(CPMobjlist)
  }
}

##################################################

Extract_CPM_ASFR_function <- function(JOY=1936:2014,obs=obs,age1=15,age2=44,parameter= c("linear","arima"),len=50,name="CoaleTrussell",file.head="Forecastlist",full_names=full_names,folder="Results/CV/Method03_CoaleTrussell1974",out="Results/CV",remove_raw=T){
  CVpredASFR  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  cohorts     <- colnames(asfr_period_to_cohort(CVpredASFR[1,1,,]))
  CVpredCASFR <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  cat("Extract results")
  joy <- JOY[1]
  for(joy in JOY){
    if(name != "")
    label <- ifelse(all(!is.na(parameter)),paste(name,paste(parameter,collapse = "_"),sep="_"),name)
    if(name == "")
    label <- ifelse(all(!is.na(parameter)),paste(paste(parameter,collapse = "_"),sep="_"),name)
    
    flag <- paste(file.head,label,"obs",obs,"joy",joy,"len",len,sep="_")
    file <- file.path(folder,paste(flag,".Rdata",sep=""))
    load(file)
    #cat(joy," ")
    for(pop in names(full_names)){
      #print(pop)
      
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[4]])
      #print(temp)
      #colnames(temp) <- paste(age1:age2)
      
      if(any(!is.na(temp))){
        CVpredASFR[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- temp[,intersect(colnames(temp),paste(1895:2050))]
      }
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[5]])
      #print(temp)
      if(any(!is.na(temp))){
        CVpredCASFR[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- temp[,intersect(colnames(temp),paste(cohorts))]
      }
    }
    
    ####
    if(remove_raw){
      unlink(file, recursive = FALSE, force = FALSE)
    }
  }
  
  flag <- paste("Forecastlist",label,"obs",obs,"len",len,sep="_")
  file <- file.path(out,paste(flag,"_CVpredASFR.Rdata",sep=""))
  save(CVpredASFR,file = file)
  file <- file.path(out,paste(flag,"_CVpredCASFR.Rdata",sep=""))
  save(CVpredCASFR,file = file)
  cat("\n")
}

##############################################################
Extract_CPM_ASFR_PI_function <- function(JOY=1936:2014,obs=obs,age1=15,age2=44,parameter= c("linear","arima"),len=50,name="CoaleTrussell",file.head="Forecastlist",full_names=full_names,folder="Results/CV/Method03_CoaleTrussell1974",out="Results/CV"){
  CVpredASFR  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  cohorts     <- colnames(asfr_period_to_cohort(CVpredASFR[1,1,,]))
  CVpredCASFR <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  
  CVpredASFRlowerPI80  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  CVpredASFRlowerPI95  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  CVpredASFRupperPI80  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  CVpredASFRupperPI95  <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  cohorts     <- colnames(asfr_period_to_cohort(CVpredASFR[1,1,,]))
  CVpredCASFRlowerPI80 <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  CVpredCASFRlowerPI95 <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  CVpredCASFRupperPI80 <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  CVpredCASFRupperPI95 <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  
  joy <- JOY[1]
  for(joy in JOY){
    if(name != "")
      label <- ifelse(all(!is.na(parameter)),paste(name,paste(parameter,collapse = "_"),sep="_"),name)
    if(name == "")
      label <- ifelse(all(!is.na(parameter)),paste(paste(parameter,collapse = "_"),sep="_"),name)
    
    flag <- paste(file.head,label,"obs",obs,"joy",joy,"len",len,sep="_")
    file <- file.path(folder,paste(flag,".Rdata",sep=""))
    load(file)
    print(joy)
    for(pop in names(full_names)){
      print(pop)
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[7]])
      if(any(!is.na(temp))){
        CVpredASFRlowerPI80[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- temp[,intersect(colnames(temp),paste(1895:2050))]
      }
      
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[11]])
      if(any(!is.na(temp))){
        CVpredASFRlowerPI95[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- temp[,intersect(colnames(temp),paste(1895:2050))]
      }
      
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[9]])
      if(any(!is.na(temp))){
        CVpredASFRupperPI80[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- temp[,intersect(colnames(temp),paste(1895:2050))]
      }
      
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[13]])
      if(any(!is.na(temp))){
        CVpredASFRupperPI95[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- temp[,intersect(colnames(temp),paste(1895:2050))]
      }
      
      ########
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[6]])
      if(any(!is.na(temp))){
        CVpredCASFRlowerPI80[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- temp[,intersect(colnames(temp),paste(cohorts))]
      }
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[10]])
      if(any(!is.na(temp))){
        CVpredCASFRlowerPI95[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- temp[,intersect(colnames(temp),paste(cohorts))]
      }
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[8]])
      if(any(!is.na(temp))){
        CVpredCASFRupperPI80[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- temp[,intersect(colnames(temp),paste(cohorts))]
      }
      temp <- as.matrix(CPMobjlist[[which(names(CPMobjlist) == pop)]][[12]])
      if(any(!is.na(temp))){
        CVpredCASFRupperPI95[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- temp[,intersect(colnames(temp),paste(cohorts))]
      }
    }
  }
  
  flag <- paste("Forecastlist",label,"obs",obs,"len",len,sep="_")
  file <- file.path(out,paste(flag,"_CVpredASFRlowerPI80.Rdata",sep=""))
  save(CVpredASFRlowerPI80,file = file)
  file <- file.path(out,paste(flag,"_CVpredASFRlowerPI95.Rdata",sep=""))
  save(CVpredASFRlowerPI95,file = file)
  file <- file.path(out,paste(flag,"_CVpredASFRupperPI80.Rdata",sep=""))
  save(CVpredASFRupperPI80,file = file)
  file <- file.path(out,paste(flag,"_CVpredASFRupperPI95.Rdata",sep=""))
  save(CVpredASFRupperPI95,file = file)
  
  file <- file.path(out,paste(flag,"_CVpredCASFRlowerPI80.Rdata",sep=""))
  save(CVpredCASFRlowerPI80,file = file)
  file <- file.path(out,paste(flag,"_CVpredCASFRlowerPI95.Rdata",sep=""))
  save(CVpredCASFRlowerPI95,file = file)
  file <- file.path(out,paste(flag,"_CVpredCASFRupperPI80.Rdata",sep=""))
  save(CVpredCASFRupperPI80,file = file)
  file <- file.path(out,paste(flag,"_CVpredCASFRupperPI95.Rdata",sep=""))
  save(CVpredCASFRupperPI95,file = file)
}
##################################################

Extract_CPM_Observed_ASFR_function <- function(asfr_period_hfd,JOY=1936:2014,obs=40,age1=15,age2=44,len=50,full_names=full_names,out="Results/"){
  CVpredASFR <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(1895:2050)),dimnames=list(names(full_names),1936:2014,age1:age2,1895:2050))
  cohorts <- colnames(asfr_period_to_cohort(CVpredASFR[1,1,,]))
  CVpredCASFR <- array(NA,dim=c(length(full_names),length(1936:2014),length(age1:age2),length(cohorts)),dimnames=list(names(full_names),1936:2014,age1:age2,cohorts))
  
  for(joy in JOY){
    print(joy)
    for(pop in names(full_names)){
      print(pop)
      ASFR <- asfr_period_hfd[paste(age1:age2),,pop]
      if(all(!is.na(ASFR[1,paste((joy-obs+1):joy)]))){
        temp <- asfr_period_hfd[paste(age1:age2),paste(JOY),pop]
        temp[,1:which(colnames(temp)==joy)] <- NA
        temp[,intersect(which(colnames(temp)==joy)+c(len:200),c(1:ncol(temp)))] <- NA
        
        CVpredASFR[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(1895:2050))] <- as.matrix(temp[,intersect(colnames(temp),paste(1895:2050))])
        temp <- asfr_period_to_cohort(temp)
        CVpredCASFR[pop,paste(joy),rownames(temp),intersect(colnames(temp),paste(cohorts))]  <- as.matrix(temp[,intersect(colnames(temp),paste(cohorts))])
      }
    }
  }
  
  flag <- paste("Forecastlist_RealData","obs",obs,"len",len,sep="_")
  file <- file.path(out,paste(flag,"_CVpredASFR.Rdata",sep=""))
  save(CVpredASFR,file = file)
  file <- file.path(out,paste(flag,"_CVpredCASFR.Rdata",sep=""))
  save(CVpredCASFR,file = file)
}
########################################################

Extract_CPM_RealPred_Given_JOY_function <- function(asfr_period_hfd,obs=5,age1=15,age2=44,parameter= c(5),len=50,name="Myrskyla2013",file.head="CVlist",folder="Results/CV/Method18_Myrskyla2013",out="Results/CV"){
  Years <- as.numeric(dimnames(asfr_period_hfd)[[2]])
  Country <- dimnames(asfr_period_hfd)[[3]]
  CPMobj <- list()
  for(pop in Country){
    joy <- max(Years[!is.na(asfr_period_hfd[1,,pop])])
    
    if(name != "")
      label <- ifelse(all(!is.na(parameter)),paste(name,paste(parameter,collapse = "_"),sep="_"),name)
    if(name == "")
      label <- ifelse(all(!is.na(parameter)),paste(paste(parameter,collapse = "_"),sep="_"),name)
    
    flag <- paste(file.head,label,"obs",obs,"joy",joy,"len",len,sep="_")
    file <- file.path(folder,paste(flag,".Rdata",sep=""))
    load(file)
    print(paste(pop,joy))
    
    CPMobj[[pop]] <- CPMobjlist[[pop]]
  }
  
  write_CPMobjlist_to_file(CPMobj,out,header="Forecastlist")
}


