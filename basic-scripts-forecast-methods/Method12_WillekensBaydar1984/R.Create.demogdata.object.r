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
#    2016.09.23
#    
#    Method: Hyndman and Ullah 2007
#     
#    Function.1: Create demogdata object from raw data matrices
#    R.Create.demogdata.object.r
############################################################
#### Input file: 
#    raw.data 
#        Matrix of data: either mortality rates or fertility rates
#### Parameters
#      year1: start year of sample data 
#      year2: end year of sample data
#      age1 : star age
#      age2 : end age
#      type : Character string showing type of demographic series: either "mortality", "fertility" or "migration".
#      label: Name of area from which the data are taken.
#      name : Name of series: usually male, female or total.
#      lambda : Box-Cox transformation parameter.
#### Output: obj.demo
#      year Vector of years
#      age Vector of ages
#      rate A list containing one or more rate matrices with one age group per row and one column per year.
#      pop A list of the same form as rate but containing population numbers instead of demographic rates.
#      type Type of object: "mortality", "fertility" or "migration".
#      label label
#      lambda lambda
#### NOTE:
#      
#### Ref:
#        
#        
###############################################################################
#install.packages("demography")
require(demography)
################################
R.Create.demogdata.object.r <- function(raw.data, 
                                        year1 = 1895, 
                                        year2 = 2014, 
                                        age1 = 15, 
                                        age2 = 44,
                                        type = "fertility",
                                        label= "test",
                                        name = "total",
                                        lambda = 0
                                        ){
  #### extract data
  age <- as.numeric(rownames(raw.data))
  period <- as.numeric(colnames(raw.data))
  
  data <- raw.data[age %in% c(age1:age2), period %in% c(year1:year2)]
  
  if(0>1){
    #### data adjust 
    if(max(data) > 50 & lambda ==0){
      CONST <- 0.001
      data <- data + CONST
    }
    if(max(data) < 20 & lambda ==0){
      CONST <- 0.000001
      data <- data + CONST
    }
  }

  #### 
  age <- as.numeric(rownames(data))
  period <- as.numeric(colnames(data))
  pop <- matrix(1000000,nrow = length(age), ncol= length(period))
  dimnames(pop) <-dimnames(data)
  
  obj.demo <- demogdata(data, pop, age, period, type, label, name, lambda)
  return(obj.demo)
}


########
if(0>1){
  pop  = "USA"
  joy  = 1975
  ASFR <- asfr_period_hfd[,,pop]
  
  raw.data <- ASFR
  year1 = joy-24 
  year2 = joy 
  age1 = 15
  age2 = 44
  VARNAME = "ASFR"
  obj.demo <- R.Create.demogdata.object.r(raw.data, year1=year1, year2=year2,age1=age1,age2=age2,lambda=0.5)
  obj.demo <- R.Create.demogdata.object.r(raw.data, year1=1895, year2=2014,age1=15,age2=44,lambda=0.5)
}
