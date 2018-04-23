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

######################################################################################################
#
#
######################################################################################################
###############################################
#### Run projection for each JOY
###############################################
run_projection_for_selected_methods <- function(asfr_period_data,expos_period_data,JOYs=1990:2014,validation_methods_info,forecasts_by_method_folder="Results/CV_ASFR_projection"){

  if(!dir.exists(forecasts_by_method_folder)){
    dir.create(forecasts_by_method_folder) 
  }
  method_names <- names(validation_methods_info)
  
  full_names <- dimnames(asfr_period_data)[[3]]
  names(full_names) <- full_names
  #####################
  #Method00_FreezeRates
  if("Method00_FreezeRates"            %in% method_names ){ 
    info <- validation_methods_info$"Method00_FreezeRates"
    
    if(info$CV==TRUE){
      cat("Method00_FreezeRates: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method00_FreezeRates")
      
      source("basic-scripts-forecast-methods/Method00_FreezeRates/Method00_FreezeRates.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method00_FreezeRates.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder,name="FreezeRates")

      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="FreezeRates")
    }
  }
  
  
  #####################
  #Method01_Hadwiger1940
  if("Method01_Hadwiger1940"            %in% method_names ){ 
    info <- validation_methods_info$"Method01_Hadwiger1940"
    
    if(info$CV==TRUE){
      cat("Method01_Hadwiger1940: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("original")
      }
      folder=file.path(forecasts_by_method_folder, "Method01_Hadwiger1940")
      
      source("basic-scripts-forecast-methods/Method01_Hadwiger1940/Method01_Hadwiger1940.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method01_Hadwiger1940.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Hadwiger")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Hadwiger")
    }
  }
  
  
  #####################
  #Method02_CoaleMcNeil1972
  if("Method02_CoaleMcNeil1972"            %in% method_names ){ 
    info <- validation_methods_info$"Method02_CoaleMcNeil1972"
    
    if(info$CV==TRUE){
      cat("Method02_CoaleMcNeil1972: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("original")
      }
      folder=file.path(forecasts_by_method_folder, "Method02_CoaleMcNeil1972")
      
      source("basic-scripts-forecast-methods/Method02_CoaleMcNeil1972/Method02_CoaleMcNeil1972.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method02_CoaleMcNeil1972.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="CoaleMcNeil")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="CoaleMcNeil")
    }
  } 
  
  
  #####################
  #Method03_CoaleTrussell1974
  if("Method03_CoaleTrussell1974"            %in% method_names ){ 
    info <- validation_methods_info$"Method03_CoaleTrussell1974"
    
    if(info$CV==TRUE){
      cat("Method03_CoaleTrussell1974: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter= c("linear","freeze")
      }
      folder=file.path(forecasts_by_method_folder, "Method03_CoaleTrussell1974")
      
      source("basic-scripts-forecast-methods/Method03_CoaleTrussell1974/Method03_CoaleTrussell1974.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method03_CoaleTrussell1974.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="CoaleTrussell")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="CoaleTrussell")
    }
  }
  
  #####################
  #Method04_Brass1974
  if("Method04_Brass1974"            %in% method_names ){ 
    info <- validation_methods_info$"Method04_Brass1974"
    
    if(info$CV==TRUE){
      cat("Method04_Brass1974: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter= c("freeze","arima","CT")
      }
      folder=file.path(forecasts_by_method_folder, "Method04_Brass1974")
      
      source("basic-scripts-forecast-methods/Method04_Brass1974/Method04_Brass1974.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method04_Brass1974.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Brass")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Brass")
    }
  }
  
  #####################
  #Method05_Evans1986
  if("Method05_Evans1986"            %in% method_names ){ 
    info <- validation_methods_info$"Method05_Evans1986"
    
    if(info$CV==TRUE){
      cat("Method05_Evans1986: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter="imputed"
      }
      folder=file.path(forecasts_by_method_folder, "Method05_Evans1986")
      
      source("basic-scripts-forecast-methods/Method05_Evans1986/Method05_Evans1986.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method05_Evans1986.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Evans")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Evans")
    }
  }
  
  #####################
  #Method06_Chandola1999
  if("Method06_Chandola1999"            %in% method_names ){ 
    info <- validation_methods_info$"Method06_Chandola1999"
    
    if(info$CV==TRUE){
      cat("Method06_Chandola1999: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("original")
      }
      folder=file.path(forecasts_by_method_folder, "Method06_Chandola1999")
      
      source("basic-scripts-forecast-methods/Method06_Chandola1999/Method06_Chandola1999.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method06_Chandola1999.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Chandola")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Chandola")
    }
  }
  
  #####################
  #Method07_Schmertmann2003
  if("Method07_Schmertmann2003"            %in% method_names ){ 
    info <- validation_methods_info$"Method07_Schmertmann2003"
    
    if(info$CV==TRUE){
      cat("Method07_Schmertmann2003: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter="freeze"
      }
      folder=file.path(forecasts_by_method_folder, "Method07_Schmertmann2003")
      
      source("basic-scripts-forecast-methods/Method07_Schmertmann2003/Method07_Schmertmann2003.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method07_Schmertmann2003.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Schmertmann2003")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Schmertmann2003")
    }
  }
  
  #####################
  #Method08_PeristeraKostaki2007M1
  if("Method08_PeristeraKostaki2007M1"            %in% method_names ){ 
    info <- validation_methods_info$"Method08_PeristeraKostaki2007M1"
    
    if(info$CV==TRUE){
      cat("Method08_PeristeraKostaki2007M1: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("original")
      }
      folder=file.path(forecasts_by_method_folder, "Method08_PeristeraKostaki2007M1")
      
      source("basic-scripts-forecast-methods/Method08_PeristeraKostaki2007M1/Method08_PeristeraKostaki2007M1.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method08_PeristeraKostaki2007M1.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="PeristeraKostakiM1")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="PeristeraKostakiM1")
    }
  }
  
  #####################
  #Method09_PeristeraKostaki2007M2
  if("Method09_PeristeraKostaki2007M2"            %in% method_names ){ 
    info <- validation_methods_info$"Method09_PeristeraKostaki2007M2"
    
    if(info$CV==TRUE){
      cat("Method09_PeristeraKostaki2007M2: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("original")
      }
      folder=file.path(forecasts_by_method_folder, "Method09_PeristeraKostaki2007M2")
      
      source("basic-scripts-forecast-methods/Method09_PeristeraKostaki2007M2/Method09_PeristeraKostaki2007M2.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method09_PeristeraKostaki2007M2.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="PeristeraKostakiM2")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="PeristeraKostakiM2")
    }
  }
  
  #####################
  #Method10_MyrskylaGoldstein2013
  if("Method10_MyrskylaGoldstein2013"            %in% method_names ){ 
    info <- validation_methods_info$"Method10_MyrskylaGoldstein2013"
    
    if(info$CV==TRUE){
      cat("Method10_MyrskylaGoldstein2013: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("arima",F,"param")
      }
      folder=file.path(forecasts_by_method_folder, "Method10_MyrskylaGoldstein2013")
      
      source("basic-scripts-forecast-methods/Method10_MyrskylaGoldstein2013/Method10_MyrskylaGoldstein2013.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method10_MyrskylaGoldstein2013.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="MG2013")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="MG2013")
    }
  }
  
  #####################
  #Method11_Saboia1977
  if("Method11_Saboia1977"            %in% method_names ){ 
    info <- validation_methods_info$"Method11_Saboia1977"
    
    if(info$CV==TRUE){
      #cat("Method11_Saboia1977: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter= c("linear","freeze")
      }
      folder=file.path(forecasts_by_method_folder, "Method11_Saboia1977")
      
      #source("basic-scripts-forecast-methods/Method11_Saboia1977/Method11_Saboia1977.R")
      for(joy in JOYs){ 
        #Forecast_CPM_ASFR_function(Method11_Saboia1977.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Saboia")

    }
  }
  
  #####################
  #Method12_WillekensBaydar1984
  if("Method12_WillekensBaydar1984"            %in% method_names ){ 
    info <- validation_methods_info$"Method12_WillekensBaydar1984"
    
    if(info$CV==TRUE){
      cat("Method12_WillekensBaydar1984: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method12_WillekensBaydar1984")
      
      source("basic-scripts-forecast-methods/Method12_WillekensBaydar1984/Method12_WillekensBaydar1984.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method12_WillekensBaydar1984.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="WillekensBaydar")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="WillekensBaydar")
    }
  }
  
  
  #####################
  #Method13_deBeer1985and1989
  if("Method13_deBeer1985and1989"            %in% method_names ){ 
    info <- validation_methods_info$"Method13_deBeer1985and1989"
    
    if(info$CV==TRUE){
      cat("Method13_deBeer1985and1989: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c(1,0,0,2,1,0)
      }
      folder=file.path(forecasts_by_method_folder, "Method13_deBeer1985and1989")
      
      source("basic-scripts-forecast-methods/Method13_deBeer1985and1989/Method13_deBeer1985and1989.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method13_deBeer1985and1989.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="deBeer")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="deBeer")
    }
  }
  
  #####################
  #Method14_Lee1993Log
  if("Method14_Lee1993Log"            %in% method_names ){ 
    info <- validation_methods_info$"Method14_Lee1993Log"
    
    if(info$CV==TRUE){
      cat("Method14_Lee1993Log: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("non","actual","arma")
      }
      folder=file.path(forecasts_by_method_folder, "Method14_Lee1993Log")
      
      source("basic-scripts-forecast-methods/Method14_Lee1993Log/Method14_Lee1993Log.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_EXPOS_function(Method14_Lee1993Log.R,asfr_period_data,expos_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="LeeLog")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="LeeLog")
    }
  }
  
  #####################
  #Method15_KohlerOrtega2002
  if("Method15_KohlerOrtega2002"            %in% method_names ){ 
    info <- validation_methods_info$"Method15_KohlerOrtega2002"
   
    if(info$CV==TRUE){
      #cat("Method15_KohlerOrtega2002: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method15_KohlerOrtega2002")
      
      #source("basic-scripts-forecast-methods/Method15_KohlerOrtega2002/Method15_KohlerOrtega2002_Load_data_from_HFD.R")
      #source("basic-scripts-forecast-methods/Method15_KohlerOrtega2002/Method15_KohlerOrtega2002_Forecast_CPM.R")
      #source("basic-scripts-forecast-methods/Method15_KohlerOrtega2002/Method15_KohlerOrtega2002.R")
      for(joy in JOYs){ 
        #Forecast_CPM_ASFR_function(Method15_KohlerOrtega2002.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="HyndmanUllah")

    }
  }

  #####################
  #Method16_HyndmanUllah2007
  if("Method16_HyndmanUllah2007"            %in% method_names ){ 
    info <- validation_methods_info$"Method16_HyndmanUllah2007"
    
    if(info$CV==TRUE){
      cat("Method16_HyndmanUllah2007: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("theoretical","3","classical","arima","noadj")
      }
      folder=file.path(forecasts_by_method_folder, "Method16_HyndmanUllah2007")
      
      source("basic-scripts-forecast-methods/Method16_HyndmanUllah2007/Method16_HyndmanUllah2007.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_EXPOS_function(Method16_HyndmanUllah2007.R,asfr_period_data,expos_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="HyndmanUllah")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="HyndmanUllah")
    }
  }
  
  #####################
  #Method17_ChengLin2010
  if("Method17_ChengLin2010"            %in% method_names ){ 
    info <- validation_methods_info$"Method17_ChengLin2010"
    
    if(info$CV==TRUE){
      cat("Method17_ChengLin2010: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("kps","quadratic")
      }
      folder=file.path(forecasts_by_method_folder, "Method17_ChengLin2010")
      
      source("basic-scripts-forecast-methods/Method17_ChengLin2010/Method17_ChengLin2010.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method17_ChengLin2010.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="ChengLin")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="ChengLin")
    }
  }
  
  #####################
  #Method18_Myrskyla2013
  if("Method18_Myrskyla2013"            %in% method_names ){ 
    info <- validation_methods_info$"Method18_Myrskyla2013"
    
    if(info$CV==TRUE){
      cat("Method18_Myrskyla2013: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c(5)
      }
      folder=file.path(forecasts_by_method_folder, "Method18_Myrskyla2013")
      
      source("basic-scripts-forecast-methods/Method18_Myrskyla2013/Method18_Myrskyla2013.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Method18_Myrskyla2013.R(asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder, file.head = "CVlist",name="Myrskyla2013")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Myrskyla2013")
    }
  }

  
  #####################
  #Method19_BayesTFR
  if("Method19_BayesTFR"            %in% method_names ){ 
    info <- validation_methods_info$"Method19_BayesTFR"
    
    if(info$CV==TRUE){
      #cat("Method19_BayesTFR: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method19_BayesTFR")
      
      #source("basic-scripts-forecast-methods/Method19_BayesTFR/Method19_BayesTFR.R")
      for(joy in JOYs){ 
        #Forecast_CPM_ASFR_function(Method19_BayesTFR.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="BayesTFR")
      
    }
  }
  
  #####################
  #Method20_BayesUN
  if("Method20_BayesUN"            %in% method_names ){ 
    info <- validation_methods_info$"Method20_BayesUN"
    
    if(info$CV==TRUE){
      #cat("Method20_BayesUN: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=c("SPE","EQU")
      }
      folder=file.path(forecasts_by_method_folder, "Method20_BayesUN")
      
      #source("basic-scripts-forecast-methods/Method20_BayesUN/Method20_BayesUN.R")
      for(joy in JOYs){ 
        #Method20_BayesUN.R(Method20_BayesUN.R,asfr_period_data,expos_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="BayesUN")

    }
  }
  

  #####################
  #Method21_Schmertmann2014
  if("Method21_Schmertmann2014"            %in% method_names ){ 
    info <- validation_methods_info$"Method21_Schmertmann2014"
    
    if(info$CV==TRUE){
      #cat("Method21_Schmertmann2014: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method21_Schmertmann2014")
      
      #source("basic-scripts-forecast-methods/Method21_Schmertmann2014/Method21_Schmertmann2014.R")
      for(joy in JOYs){ 
        #Forecast_CPM_ASFR_function(Method21_Schmertmann2014.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #####
    }
  }
  
  #####################
  #Method22_LiWu2003
  if("Method22_LiWu2003"            %in% method_names ){ 
    info <- validation_methods_info$"Method22_LiWu2003"
    
    if(info$CV==TRUE){
      cat("Method22_LiWu2003: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method22_LiWu2003")
      
      source("basic-scripts-forecast-methods/Method22_LiWu2003/Method22_LiWu2003.R")
      for(joy in JOYs){ 
        cat(joy," ")
        Forecast_CPM_ASFR_function(Method22_LiWu2003.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="LiWu")
      save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="LiWu")
    }
  }

  #####################
  #Method23_Sobotka2011
  if("Method23_Sobotka2011"            %in% method_names ){ 
    info <- validation_methods_info$"Method23_Sobotka2011"
    
    if(info$CV==TRUE){
      #cat("Method23_Sobotka2011: ")
      obs <- info$obs
      age1 <- info$ages[1]
      age2 <- info$ages[2]
      len <- info$len
      if(info$parameter!="default"){
        parameter= info$parameter
      }
      if(info$parameter=="default"){
        parameter=NA
      }
      folder=file.path(forecasts_by_method_folder, "Method23_Sobotka2011")
      
      #source("basic-scripts-forecast-methods/Method23_Sobotka2011/Method23_Sobotka2011.R")
      for(joy in JOYs){ 
        #Forecast_CPM_ASFR_function(Method23_Sobotka2011.R,asfr_period_data,joy=joy,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,out=F)
      }
      
      #Extract_CPM_ASFR_function(JOY=JOYs,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,full_names=full_names,out = folder ,name="Sobotka")
      #save_observed_predicted_cfr(asfr_period_data,obs=obs,age1=age1,age2=age2,parameter=parameter,len=len,folder=folder,name="Sobotka")
    }
  }

  ##############
}
