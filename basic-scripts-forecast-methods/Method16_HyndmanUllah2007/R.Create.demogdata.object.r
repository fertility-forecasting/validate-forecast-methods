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
#    rate.data 
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
#require(demography)
####
demogdata <- function (data, pop, ages, years, type, label, name, lambda) 
{
  p <- nrow(data)
  n <- ncol(data)
  if (nrow(pop) != p | ncol(pop) != n) 
    stop("data and pop are of different size")
  if (length(ages) != p) 
    stop("Number of ages doesn't match data")
  if (length(years) != n) 
    stop("Number of years doesn't match data")
  types <- c("mortality", "fertility", "migration", "population")
  idx <- pmatch(type, types)
  if (is.na(idx)) 
    warning("Unknown type")
  else type <- types[idx]
  if (missing(lambda)) {
    if (type == "mortality") 
      lambda <- 0
    else if (type == "fertility") 
      lambda <- 0.4
    else lambda <- 1
  }
  if (type == "population") {
    obj <- list(year = years, age = ages, pop = list(as.matrix(pop)), 
                type = type, label = label, lambda = lambda)
    dimnames(obj$pop[[1]]) <- list(ages, years)
  }
  else {
    obj <- list(year = years, age = ages, rate = list(as.matrix(data)), 
                pop = list(as.matrix(pop)), type = type, label = label, 
                lambda = lambda)
    dimnames(obj$rate[[1]]) <- dimnames(obj$pop[[1]]) <- list(ages, 
                                                              years)
    names(obj$rate) <- names(obj$pop) <- name
  }
  return(structure(obj, class = "demogdata"))
}

####
extract.years <- function (data, years) 
{
  idx <- match(years, data$year)
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) 
    stop("No data available for those years")
  no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if (!no.rate) {
    nn <- length(data$rate)
    dname <- dimnames(data$rate[[1]])
    dname <- list(dname[[1]], dname[[2]][idx])
  }
  else if (!no.pop) {
    nn <- length(data$pop)
    dname <- dimnames(data$pop[[1]])
    dname <- list(dname[[1]], dname[[2]][idx])
  }
  else stop("No data!")
  for (j in 1:nn) {
    if (!no.rate) {
      data$rate[[j]] <- matrix(data$rate[[j]][, idx], ncol = length(idx))
      dimnames(data$rate[[j]]) <- dname
    }
    if (!no.pop) {
      data$pop[[j]] <- matrix(data$pop[[j]][, idx], ncol = length(idx))
      dimnames(data$pop[[j]]) <- dname
    }
    if (!is.null(data$obs.var)) {
      if (length(dim(data$obs.var[[j]])) == 2) 
        data$obs.var[[j]] <- data$obs.var[[j]][, idx]
      else data$obs.var[[j]] <- data$obs.var[[j]][idx]
    }
  }
  data$year <- data$year[idx]
  return(data)
}

####
extract.ages <- function (data, ages, combine.upper = TRUE) 
{
  no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if (combine.upper) {
    if (no.pop) {
      if (max(data$age) > max(ages)) 
        warning("No population data available for combining upper ages")
    }
    else data <- set.upperage(data, max(ages))
  }
  idx <- match(ages, data$age)
  idx <- idx[!is.na(idx)]
  no.pop <- is.null(data$pop)
  no.rate <- is.null(data$rate)
  if (!no.rate) {
    nn <- length(data$rate)
    dname <- dimnames(data$rate[[1]])
    dname <- list(dname[[1]][idx], dname[[2]])
  }
  else if (!no.pop) {
    nn <- length(data$pop)
    dname <- dimnames(data$pop[[1]])
    dname <- list(dname[[1]][idx], dname[[2]])
  }
  else stop("No data!")
  for (j in 1:nn) {
    if (!no.rate) {
      data$rate[[j]] <- matrix(data$rate[[j]][idx, ], nrow = length(idx))
      dimnames(data$rate[[j]]) <- dname
    }
    if (!no.pop) {
      data$pop[[j]] <- matrix(data$pop[[j]][idx, ], nrow = length(idx))
      dimnames(data$pop[[j]]) <- dname
    }
    if (!is.null(data$obs.var)) {
      if (length(dim(data$obs.var[[j]])) == 2) 
        data$obs.var[[j]] <- data$obs.var[[j]][idx, ]
      else data$obs.var[[j]] <- data$obs.var[[j]][idx]
    }
  }
  data$age <- data$age[idx]
  return(data)
}

################################
R.Create.demogdata.object.r <- function(rate.data,
                                        pop.data,
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
  age <- as.numeric(rownames(rate.data))
  period <- as.numeric(colnames(rate.data))
  data <- rate.data[age %in% c(age1:age2), period %in% c(year1:year2)]
  
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
  age <- as.numeric(rownames(pop.data))
  period <- as.numeric(colnames(pop.data))
  pop <- pop.data[age %in% c(age1:age2), period %in% c(year1:year2)]

  ####
  age <- as.numeric(rownames(data))
  period <- as.numeric(colnames(data))
  obj.demo <- demogdata(data, pop, age, period, type, label, name, lambda)
  return(obj.demo)
}


########
if(0>1){
  pop  = "USA"
  joy  = 1975
  ASFR <- asfr_period_hfd[,,pop]
  
  rate.data <- ASFR
  year1 = joy-24 
  year2 = joy 
  age1 = 15
  age2 = 44
  VARNAME = "ASFR"
  obj.demo <- R.Create.demogdata.object.r(rate.data, year1=year1, year2=year2,age1=age1,age2=age2,lambda=0.5)
  obj.demo <- R.Create.demogdata.object.r(rate.data, year1=1895, year2=2014,age1=15,age2=44,lambda=0.5)
}
