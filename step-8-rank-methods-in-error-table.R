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
#################################################################

error.summary <- function(x){
  RankStat <- c(quantile(x,prob=c(0.50,0.80,0.90,0.95),na.rm=TRUE),mean(x,na.rm=T),sqrt(mean(x^2,na.rm=T)))
  return(round(RankStat,4))
}

#### function to calculat Counting Inversion of ranks comparing to our ranking as reference.

CountingInversion <- function(x){
  n <- length(x)
  ct <- 0
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      if(x[i]>x[j]){
        ct <- ct +1
        #print(sort(c(x[i],x[j])))
      }
    }
  }
  return(ct)
}
#ct <- CountingInversion(c(2,4,1,3,5))
##################################
APE_ks_test <- function(obj.list){
  obj.list.ks <- diag(length(obj.list))
  rownames(obj.list.ks) <- substr(names(obj.list),10,33)
  colnames(obj.list.ks) <- substr(names(obj.list),10,33)
  for(i in 1:length(obj.list)){
    for(j in 1:length(obj.list)){
      #print(c(i,j))
      ks.obj <- ks.test(obj.list[[i]], obj.list[[j]],alternative = "greater")
      obj.list.ks[i,j] <- ks.obj$p.value
    }
  }
  return(obj.list.ks)
}


############################################ 
CV_APE_RankStat_table <- function(input,file="RankStat_table.txt"){
  RankStat <- NULL

  
  for(i in 1:length(input)){
    
    RankStat <- rbind(RankStat, error.summary(input[[i]]*100))
  }
  
  ks_obj <- APE_ks_test(input)
  
  ks <- ks_obj < 0.05 & t(ks_obj >= 0.05)
  ks.rank <- apply(ks,1,sum)
  RankStat <- cbind(RankStat, ks.rank)
  RK <- rank(RankStat[,5]) + 10*rank(-RankStat[,7])
  
  RankStat <- cbind(RankStat,RK)
  colnames(RankStat)[5:8] <- c("Mean","RMSE","KS","Rank")
  RankStat
  
  tab <- cbind(paste(rank(RK,ties.method = "min"),rownames(RankStat),sep=". "),
               paste(round(RankStat[,7],0)," (",rank(-RankStat[,7],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,1],2)," (",rank(RankStat[,1],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,2],2)," (",rank(RankStat[,2],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,3],2)," (",rank(RankStat[,3],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,4],2)," (",rank(RankStat[,4],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,5],2)," (",rank(RankStat[,5],ties.method = "min"),")",sep=""),
               paste(round(RankStat[,6],2)," (",rank(RankStat[,6],ties.method = "min"),")",sep=""))
  
  require(xtable)
  colnames(tab) <- c("Method","KS","50%","80%","90%","95%","Mean","RMSE")
  rownames(tab) <- names(input)
  
  ############################
  RankAll <- cbind(rank(RK,ties.method = "min"),
                   rank(-RankStat[,7],ties.method = "min"),
                   rank(RankStat[,1],ties.method = "min"),
                   rank(RankStat[,2],ties.method = "min"),
                   rank(RankStat[,3],ties.method = "min"),
                   rank(RankStat[,4],ties.method = "min"),
                   rank(RankStat[,5],ties.method = "min"),
                   rank(RankStat[,6],ties.method = "min"))
  RankAll <- RankAll[order(RK),]
  
  
  Inversion_number=apply(RankAll,2,CountingInversion)
  Inversion_number[1] <- "Inversion_number"
  error_table <- rbind(tab[order(RK),],Inversion_number)
  
  print(xtable(error_table), include.rownames=FALSE,file=file)
  
}





