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



##########################################################################
#### Step 0. Set up working directory
rm(list=ls())

project.directory <- c(".")
setwd(project.directory)
##########################################################################
#### Step 1. Install and load required R packages
source("step-1-load-required-R-packages.R")

##########################################################################
#### Step 2.  Load raw data that is required to derive fertility forecasts.
####          Note that you need to provide input data yourself in the folder "input-data".
####          The input data should be in the format: age times calendar year. 
####          Note that some forecast methods require both fertility rates and exposure data. 
####          Note that the Human Fertility Database (HFD) provides fertility data (e.g. asfrRR, exposRR) 
####          at https://www.humanfertility.org/

source("step-2-load-input-data.R")

## ASFR: Age-specific fertility rates
asfr_period <- read.hfd("input-data/asfrRR.txt",vASFR = c("ASFR"),ages=c(15:49))

## Exposure: Female population exposure
expos_period <- read.hfd("input-data/exposRR.txt",vASFR = c("Exposure"),ages=c(15:49))

## Name of the regions
#full_names <- dimnames(asfr_period)[3]
#names(full_names) <- full_names
##########################################################################
#### Step 3.  Select a forecast method (as shown in file CPM_Github_Methods_Info.xlsx)
####          Please go to the following script to (1) select the forecast methods that enter the cross validation
####          and to (2) set / specify their parameters.
####
source("step-3-select-forecast-methods-and-specify-their-parameters.R")

##########################################################################
#### Step 4.  Derive fertility forecasts with all the selected forecast methods (step 3) 
####          with the function "run_projection_for_selected_methods".
####          You need to specify the input data, jump-off years, and the folder to store the forecast results.
####          Before you can run "run_projection_for_selected_methods", please run scripts 4a, 4b, and 4c as pre-steps. 
#### Step 4a. Basic functions for saving and reading forecast results.           
#### Step 4b. Run projection for selected methods.

#### 
source("step-4a-function-CV-and-Forecast-CPM.R")
source("step-4b-run-projection.R")

#### Set folde to save results

run_projection_for_selected_methods(asfr_period,
                                    expos_period,
                                    JOYs=2013:2014,
                                    validation_methods_info,
                                    forecasts_by_method_folder="results/forecasts-by-method-test")
  
##########################################################################
#### Step 5. Calculate forecast errors, absolute percentage error (APE), for selected forecast methods.
#### 
source("step-5-calculate-forecast-errors.R")

#### Please select a subset of methods for which you would like to caclulate APEs 
#### (if you have not already done so in step 3)
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
  "Method12_WillekensBaydar1984"   = T,
  "Method13_deBeer1985and1989"     = T,
  "Method14_Lee1993Log"            = T,
  "Method18_Myrskyla2013"          = T
)

#### Caclulate APEs with function "run_projection_for_selected_methods"
####

APE_list <- compute_APE_for_selected_methods(methods_for_APE,forecasts_by_method_folder="results/forecasts-by-method-test")

##########################################################################
#### Step 6. KS test and grid plot
#### 
source("step-6-visualize-statistical-test-results.R")

ks <- APE_ks_test(APE_list)
pdf(file.path(forecasts_by_method_folder,"APE_KS_test_grid_plot.pdf"),width = 8,height = 8)
APE_ks_test_grid_plot(ks)
dev.off()


##########################################################################
#### Step 7. Visualize thresholds for absolute percentage errors (APE) across forecast methods
#### 

source("step-7-visualize-thresholds-for-APE.R")
CV_APE_threshold_plot(APE_list,method.name = substr(names(APE_list),10,33),plot.file = file.path(forecasts_by_method_folder,"APE_threshold_plot.pdf"))

##########################################################################
#### Step 8. Compare and rank forecast methods based on seven error metrics   
####          (i.e. KS test statistic, 50%, 80%, 90%, 95% threshold APE, MAPE, RMSE) 

source("step-8-rank-methods-in-error-table.R")
CV_APE_RankStat_table(APE_list,file = file.path(forecasts_by_method_folder,"rank-table.txt"))




