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

###########################################################################
#### read data
read.hfd <- function(file,vYEAR="Year",vAGE="Age",vASFR="ASFR",ages=c(15:49)){
  #########################################################################
  ##### Read Raw Data
  #####################################
  #### Read from hfd: ASFR
  hfd_download_data <- list()
  
  FF <- read.table(file,skip=2,as.is=T, header=T)
  FF <- subset(FF, FF[,vAGE] %in% ages)
  
  for (k in unique(FF$Code)) {
    sel = ( (FF$Code==k)  )
    hfd_download_data[[k]] <- tapply(as.numeric(FF[sel,vASFR]), list(FF[sel,vAGE], FF[sel,vYEAR]),sum) 
  }
  
  hfd_download_data[["CHL"]] <- NULL
  #names(hfd_download_data)
  
  #########################################
  #Put ASFR in format 15-44 and 1847-2014:
  #########################################
  #print(names(hfd_download_data))
  asfr_period_hfd <- array(NA,dim=c(length(ages),length(1847:2014),length(hfd_download_data)),dimnames=list(ages,c(1847:2014),names(hfd_download_data)))
  
  for(i in 1:length(hfd_download_data)){
    #print(i)
    years <- intersect(colnames(hfd_download_data[[i]]),as.character(c(1847:2014)))
    asfr_period_hfd[as.character(ages),years,i] <- hfd_download_data[[i]][as.character(ages),years]
  }
  
  dimnames(asfr_period_hfd)[[3]] <- names(hfd_download_data)
  
  #asfr_period_hfd[,,"AUT"]
  #setwd(the.results.path)
  #save(asfr_period_hfd,file="data.asfr_period_hfd.Rdata")
  #load("data.asfr_period_hfd.Rdata")
  return(asfr_period_hfd)
}


##########################################################
asfr_period_to_cohort <- function(input){
  temp <- expand.grid(as.numeric(dimnames(input)[[1]]),as.numeric(dimnames(input)[[2]]))
  cohort <- temp[,2] - temp[,1]
  output <- tapply(as.vector(input), list(temp[,1], cohort),sum)
  return(as.data.frame(output))
}
