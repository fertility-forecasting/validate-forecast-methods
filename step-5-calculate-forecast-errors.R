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

compute_APE_for_selected_methods <- function(methods_for_APE,forecasts_by_method_folder="Results/CV_ASFR_projection"){
  
  
  methods_for_APE <- methods_for_APE[methods_for_APE]
  
  method_names <- names(methods_for_APE)
  obj.list <- as.list(rep(NA,length(method_names)))
  names(obj.list) <- method_names
  
  ind <- 1
  for(i in 1:length(obj.list)){
    name <- method_names[i]
    
    
    folder <- file.path(forecasts_by_method_folder, name)
    file_cfr_observed <- file.path(folder,"Result_cfr_CVobse.Rdata")
    file_cfr_forecast <- file.path(folder,"Result_cfr_CVpred.Rdata")
    load(file_cfr_observed)
    load(file_cfr_forecast)
    
    APE <- abs(cfr_observed - cfr_forecast)/cfr_observed
    obj.list[[i]] <- APE
    ind <- ind * APE
    
  }
  
  ind[!is.na(ind)] <- 1
  
  for(i in 1:length(obj.list)){
    obj.list[[i]] <- obj.list[[i]] * ind 
  }
  
  return(obj.list)
} 
