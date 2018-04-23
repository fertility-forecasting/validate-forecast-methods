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


########################################################################
#### "validation_methods_info" provides information about which methods enter the cross validation and their 
####  parameter settings. Please adjust the information given in this list according to your needs.
#### 1. CV=[T/F]: select or does not select a certain forecast method to be used in the cross validation
#### 2. obs:      length of observation window  
#### 3. ages:     range of ages
#### 4. len:      length of projection (/ forecast horizon)
#### 5. parameter:parameter settings used for each method. To use another parameter setting than "default" please go to...; 

validation_methods_info <- list(
  "Method00_FreezeRates"           = list(CV=T, obs=30, ages=c(15,40), len=50, parameter="default"),
  "Method01_Hadwiger1940"          = list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method02_CoaleMcNeil1972"       = list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method03_CoaleTrussell1974"     = list(CV=T, obs=40, ages=c(15,49), len=50, parameter="default"),
  "Method04_Brass1974"             = list(CV=T, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method05_Evans1986"             = list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method06_Chandola1999"          = list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method07_Schmertmann2003"       = list(CV=T, obs=40, ages=c(15,40), len=50, parameter="default"),
  "Method08_PeristeraKostaki2007M1"= list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method09_PeristeraKostaki2007M2"= list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method10_MyrskylaGoldstein2013" = list(CV=T, obs=35, ages=c(15,45), len=50, parameter="default"),
  "Method11_Saboia1977"            = list(CV=F, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method12_WillekensBaydar1984"   = list(CV=T, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method13_deBeer1985and1989"     = list(CV=T, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method14_Lee1993Log"            = list(CV=T, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method15_KohlerOrtega2002"      = list(CV=F, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method16_HyndmanUllah2007"      = list(CV=T, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method17_ChengLin2010"          = list(CV=T, obs=30, ages=c(15,44), len=50, parameter="default"),
  "Method18_Myrskyla2013"          = list(CV=T, obs=30, ages=c(15,40), len=50, parameter="default"),
  "Method19_BayesTFR"              = list(CV=F, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method20_BayesUN"               = list(CV=F, obs=30, ages=c(15,40), len=50, parameter="default"),
  "Method21_Schmertmann2014"       = list(CV=F, obs=40, ages=c(15,44), len=50, parameter="default"),
  "Method22_LiWu2003"              = list(CV=T, obs=30, ages=c(15,40), len=50, parameter="default"),
  "Method23_Sobotka2011"           = list(CV=F, obs=40, ages=c(15,44), len=50, parameter="default"))

########################################################################
#### "method_names" contains the names of the methods that enter the cross validation 
####

method_names <- names(validation_methods_info)

########################################################################
#### "methods_for_APE" selects the methods for which forecast errors will be calculated. 
####

methods_for_APE <- c(
  "Method00_FreezeRates"           = T,
  "Method01_Hadwiger1940"          = T,
  "Method02_CoaleMcNeil1972"       = T,
  "Method03_CoaleTrussell1974"     = T,
  "Method04_Brass1974"             = T,
  "Method05_Evans1986"             = T,
  "Method06_Chandola1999"          = T,
  "Method07_Schmertmann2003"       = T,
  "Method08_PeristeraKostaki2007M1"= T,
  "Method09_PeristeraKostaki2007M2"= T,
  "Method10_MyrskylaGoldstein2013" = T,
  "Method11_Saboia1977"            = F,
  "Method12_WillekensBaydar1984"   = T,
  "Method13_deBeer1985and1989"     = T,
  "Method14_Lee1993Log"            = T,
  "Method15_KohlerOrtega2002"      = F,
  "Method16_HyndmanUllah2007"      = T,
  "Method17_ChengLin2010"          = T,
  "Method18_Myrskyla2013"          = T,
  "Method19_BayesTFR"              = F,
  "Method20_BayesUN"               = F,
  "Method21_Schmertmann2014"       = F,
  "Method22_LiWu2003"              = T,
  "Method23_Sobotka2011"           = F 
)

methods_for_APE <- c(
  "Method00_FreezeRates"           = T,
  "Method01_Hadwiger1940"          = T,
  "Method02_CoaleMcNeil1972"       = T,
  "Method03_CoaleTrussell1974"     = T,
  "Method04_Brass1974"             = T,
  "Method05_Evans1986"             = T,
  "Method06_Chandola1999"          = T,
  "Method07_Schmertmann2003"       = T,
  "Method08_PeristeraKostaki2007M1"= T,
  "Method09_PeristeraKostaki2007M2"= T,
  "Method10_MyrskylaGoldstein2013" = T,
  "Method12_WillekensBaydar1984"   = T,
  "Method13_deBeer1985and1989"     = T,
  "Method14_Lee1993Log"            = T,
  "Method18_Myrskyla2013"          = T,
  "Method17_ChengLin2010"          = T,
  "Method22_LiWu2003"              = T
)
