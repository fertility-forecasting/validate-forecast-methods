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
#    2016.09.07
#    
#    Method: Cheng and Lin 2010
#     
#    Function.2: Kohler and Philipov 2001 + kernel smoothing
#    R.ChengLin2010.function2.KP.smooth.r
############################################################
#### Input file: obj.APC 
#    [output of APC framework (R.ChengLin2010.function1.APC.r) or raw data]
#    Adjusted ASFR (named as "DEPEND") will be used
#### Parameters
#    mab: mean age at birth: fixed or missing (using mean mab of each year estimated in the function) 
#### Output: obj.KPS
#      obj.KPs$CTFR.all : all CTFR estimations
#      obj.KPs$CTFR.Raw : TFR
#      obj.KPs$CTFR.bf  : TFR with BF adjustment
#      obj.KPs$CTFR.bfs : CTFR.bf + kernel smoothing
#      obj.KPs$CTFR.KP  : TFR with var adjustment (KP)
#      obj.KPs$CTFR.KPS : CTFR.KP + kernel smoothing
#      obj.KPs$mab      : mab from input or estimation (mean mab of each year)
###############################################################################
#############################################################
#### Inner function 1: smooth.43RSR2H
#     In S-Plus function smooth() uses 4(3RSR)2H method, which is
#     not included in R function smooth(). Here we write 4(3RSR)2H 
#     algorithm in R.
#### Input: series of data
#############################################################
smooth.43RSR2H <- function(x){
  x <- as.numeric(x)
  n <- length(x)
  # 4
  out = NULL
  for (i in 1:(n-1)) {out[i+1] = median(x[i:(i+3)])}
  out[is.na(out)] = x[is.na(out)]
  
  # (3RSR)
  x2 = smooth(out, "3RSR", twiceit = T)
  # 2
  out2 = NULL
  for (i in 1:(n-1)) {out2[i+1] = median(x2[i:(i+1)])}
  out2[is.na(out2)] = x2[is.na(out2)]
  # H
  out <- filter(out2, filter = c(1/4, 1/2, 1/4), sides = 2)
  out[is.na(out)] = x[is.na(out)]
  out
} 

if(0>1){
  # Data
  set.seed(1)
  x <- c(4, 1, 3, 6, 6, 4, 1, 6, 2, 4, 2)
  x <- sin(seq(1,10,0.5)) + rnorm(length(seq(1,10,0.5)),0,0.2)
  plot(x)
  lines(x)
  lines(smooth.43RSR2H(x),col="red")
  lines(smooth(x),col="green")
  lines(smooth.spline(x),col="blue")
}

#############################################################
#### Inner function 2: F.bf.calc
#    This function comes from Kohler and Philipov 2001
#    http://www.ssc.upenn.edu/~hpkohler/data-and-programs/bfvariance/bfvarianceprograms.html
#    F.bf.calc is the primary function to calculate the adjustment
#    of the parity-specific TFR with and without variance effects
#### Input
#    nfx: a matrix of ASFR, row for age and column for year
#    age: same as rows of nfx 
#    years: same as columns of nfx
#### Output
#    
#############################################################
F.bf.calc <- function(nfx, age, years, i.prec=0.001){ 
  
  # calculation of the usual TFR,
  # mean age, variance, and
  # centralized third moment
  # of the fertility schedule
  p.age <- age+0.5
  tfr <- apply(nfx, 2, sum)
  mab <- apply(nfx * (age + 0.5), 2, sum) / tfr
  vab <- apply(nfx * (age + 0.5)^2, 2, sum) / tfr - mab^2
  age.0 <- sweep(matrix(p.age,length(p.age),length(years),byrow=F),2,mab,"-")
  kub <- apply(nfx*(age.0^3), 2, sum) / tfr
  names(tfr) <- paste("yr",years)
  names(mab) <- paste("yr",years)
  names(vab) <- paste("yr",years)
  
  ## smooth the time series for TFR, mab, vab,kub
  if(0>1){
    tfr.s <- as.vector(smooth(tfr), twice=T)
    names(tfr.s) <- paste("yr",years)
    mab.s <- as.vector(smooth(mab), twice=T)
    names(mab.s) <- paste("yr",years)
    vab.s <- as.vector(smooth(vab), twice=T)
    names(vab.s) <- paste("yr",years)
    kub <- as.vector(smooth(kub), twice=T)
    names(kub) <- paste("yr",years)
  }
  
  ## using R smooth function, defult is 3RS3R
  if(0>1){
    tfr.s <- as.vector(smooth(tfr, kind="3RS3R", twiceit=T))
    names(tfr.s) <- paste("yr",years)
    mab.s <- as.vector(smooth(mab, kind="3RS3R", twiceit=T))
    names(mab.s) <- paste("yr",years)
    vab.s <- as.vector(smooth(vab, kind="3RS3R", twiceit=T))
    names(vab.s) <- paste("yr",years)
    kub <- as.vector(smooth(kub, kind="3RS3R", twiceit=T)) 
    names(kub) <- paste("yr",years)
  }
  
  ## using 4(3RSR2)H method, same as S-plus
  if(0<1){
    tfr.s <- as.vector(smooth.43RSR2H(tfr)) 
    names(tfr.s) <- paste("yr",years)
    mab.s <- as.vector(smooth.43RSR2H(mab)) 
    names(mab.s) <- paste("yr",years)
    vab.s <- as.vector(smooth.43RSR2H(vab)) 
    names(vab.s) <- paste("yr",years)
    kub <- as.vector(smooth.43RSR2H(kub))  
    names(kub) <- paste("yr",years)
  }
  # implementation of BF formula
  r.bf <- c(NA, 0.5 * (mab.s[3:length(years)] - 
                         mab.s[1:(length(years)-2)]), NA)
  tfr.adj.bf <- tfr.s / (1 - r.bf)
  names(r.bf) <- paste("yr",years)
  names(tfr.adj.bf) <- paste("yr",years)
  
  
  # calclulate delta from observed variance
  # (using the second and second to last observation of delta
  # for the first and last year of the data to avoid a reduction
  # of length of the time series)
  
  delta <- 0.25 * log(vab.s[3:length(years)]/
                        vab.s[1:(length(years)-2)])
  delta <- c(delta[1], delta, delta[length(delta)])
  names(delta) <- paste("yr",years)
  
  # calculate bias corrected bf-r
  # based on a one-step correction (Result 11)
  delta.change <- c(NA, 0.5 * (delta[3:length(years)] - 
                                 delta[1:(length(years)-2)]), NA)
  r.bf.biascorr <- r.bf + vab.s / (1 - r.bf) * (2*delta^2 + delta.change)
  
  #iterate estimation for the calculation of variance effects (Result 13)
  r.res <- F.iter.calc(mab.s,vab.s,kub,years,delta,i.prec)
  
  delta <- c(NA,delta[-c(1,length(delta))],NA)
  gamma <- r.res$gamma
  s2 <- r.res$s2
  a.bar <- r.res$a.bar
  
  # calculate adjusted TFR based on "correct" tempo gamma
  tfr.adj.var <- tfr.s/(1 - gamma)
  
  return(
    list(
    years = years,
    tfr = tfr,
    tfr.s = tfr.s,
    mab = mab,
    mab.s = mab.s,
    vab = vab,
    vab.s = vab.s,
    kub = kub,
    r.bf = r.bf,
    tfr.adj.bf = tfr.adj.bf,
    r.bf.biascorr = r.bf.biascorr,
    gamma = gamma,
    delta = delta,
    mab.adj.var = a.bar,
    vab.adj.var = s2,
    tfr.adj.var = tfr.adj.var
  )
  )
}
# F.iter.calc performs the iteration in Result 13
# and iteratively calculates gamma and then corrects
# for the distortions caused by variance effects
F.iter.calc <- function(mab,vab,kub,yrs, delta.hat, prec=0.0001,PRINT=F){
  a.hat <- mab
  s2.hat <- vab
  i <- 0
  g.diff <- 1
  
  #step 1
  gamma.hat <- 0.5*(mab[3:length(yrs)]-mab[1:(length(yrs)-2)])
  gamma.hat <- c(gamma.hat[1], gamma.hat, gamma.hat[length(gamma.hat)])
  names(gamma.hat) <- paste("yr",yrs)
  
  #iterate for calculation of gamma
  while (g.diff > prec){
    i <- i + 1
    
    #step 2
    s2.hat.old <- s2.hat
    s2.hat <- vab  +  (delta.hat/(1-gamma.hat) * s2.hat)^2 +
      delta.hat/(1-gamma.hat) * kub
    
    #step 3
    a.hat.old <- a.hat
    a.hat <- mab + delta.hat/(1-gamma.hat)*s2.hat
    
    #step 4
    gamma.hat.old <- gamma.hat
    gamma.hat <- 0.5*(a.hat[3:length(yrs)]-a.hat[1:(length(yrs)-2)])
    gamma.hat <- c(gamma.hat[1], gamma.hat, gamma.hat[length(gamma.hat)])
    # avoid estimates that are implausible
    gamma.hat[gamma.hat > 0.5] <- 0.5
    gamma.hat[gamma.hat < -0.5] <- -0.5
    names(gamma.hat) <- paste("yr",yrs)
    
    # calculate convergence criterion
    diff.vec <- abs(c(a.hat,s2.hat,gamma.hat) - 
                      c(a.hat.old,s2.hat.old,gamma.hat.old))
    g.diff <- max(diff.vec[!is.na(diff.vec)])
    
    if(PRINT){
      print(paste("iteration",i,"parameters for gamma"))
      print(g.diff)
      print(gamma.hat, digits=4)
    }

  }
  
  
  # put back missing values
  gamma.hat <- c(NA,gamma.hat[-c(1,length(gamma.hat))],NA)
  a.hat <- c(NA,a.hat[-c(1,length(a.hat))],NA)
  s2.hat <- c(NA,s2.hat[-c(1,length(s2.hat))],NA)
  
  
  names(a.hat) <- paste("yr",yrs)
  names(s2.hat) <- paste("yr",yrs)
  names(gamma.hat) <- paste("yr",yrs)
  
  list(a.bar = a.hat, s2 = s2.hat,
       gamma = gamma.hat, delta = delta.hat)
}


F.allorder <- function(nfx.tot, age, os.result){
  #direct calculation from nFx with all orders
  tfr <- apply(nfx.tot, 2, sum)
  mab <- apply(nfx.tot * (age + 0.5), 2, sum) / tfr
  vab <- apply(nfx.tot * (age + 0.5)^2, 2, sum) / tfr - mab^2
  f.years <- os.result[[1]]$years
  
  #indirect calculation from order specific results
  i.tfr <- (os.result[[1]]$tfr.s + os.result[[2]]$tfr.s +
              os.result[[3]]$tfr.s + os.result[[4]]$tfr.s)
  i.mab <- (os.result[[1]]$tfr.s * os.result[[1]]$mab.s +
              os.result[[2]]$tfr.s * os.result[[2]]$mab.s +
              os.result[[3]]$tfr.s * os.result[[3]]$mab.s +
              os.result[[4]]$tfr.s * os.result[[4]]$mab.s) / i.tfr
  i.vab <- (os.result[[1]]$tfr.s * os.result[[1]]$vab.s +
              os.result[[2]]$tfr.s * os.result[[2]]$vab.s +
              os.result[[3]]$tfr.s * os.result[[3]]$vab.s +
              os.result[[4]]$tfr.s * os.result[[4]]$vab.s) / i.tfr + 
    (os.result[[1]]$tfr.s * (os.result[[1]]$mab.s-i.mab)^2 +
       os.result[[2]]$tfr.s * (os.result[[2]]$mab.s-i.mab)^2 +
       os.result[[3]]$tfr.s * (os.result[[3]]$mab.s-i.mab)^2 +
       os.result[[4]]$tfr.s * (os.result[[4]]$mab.s-i.mab)^2) / i.tfr
  
  tfr.adj.var <- (os.result[[1]]$tfr.adj.var + 
                    os.result[[2]]$tfr.adj.var +
                    os.result[[3]]$tfr.adj.var + 
                    os.result[[4]]$tfr.adj.var)
  
  mab.adj.var <- (os.result[[1]]$tfr.s * os.result[[1]]$mab.adj.var +
                    os.result[[2]]$tfr.s * os.result[[2]]$mab.adj.var +
                    os.result[[3]]$tfr.s * os.result[[3]]$mab.adj.var +
                    os.result[[4]]$tfr.s * os.result[[4]]$mab.adj.var) / i.tfr
  
  gamma <- (os.result[[1]]$tfr.s * os.result[[1]]$gamma +
              os.result[[2]]$tfr.s * os.result[[2]]$gamma +
              os.result[[3]]$tfr.s * os.result[[3]]$gamma +
              os.result[[4]]$tfr.s * os.result[[4]]$gamma) / i.tfr
  
  delta <- (os.result[[1]]$tfr.s * os.result[[1]]$delta +
              os.result[[2]]$tfr.s * os.result[[2]]$delta +
              os.result[[3]]$tfr.s * os.result[[3]]$delta +
              os.result[[4]]$tfr.s * os.result[[4]]$delta) / i.tfr
  
  r.bf <- (os.result[[1]]$tfr.s * os.result[[1]]$r.bf +
             os.result[[2]]$tfr.s * os.result[[2]]$r.bf +
             os.result[[3]]$tfr.s * os.result[[3]]$r.bf +
             os.result[[4]]$tfr.s * os.result[[4]]$r.bf) / i.tfr
  
  vab.adj.var <- (os.result[[1]]$tfr.s * os.result[[1]]$vab.adj.var +
                    os.result[[2]]$tfr.s * os.result[[2]]$vab.adj.var +
                    os.result[[3]]$tfr.s * os.result[[3]]$vab.adj.var +
                    os.result[[4]]$tfr.s * os.result[[4]]$vab.adj.var) / i.tfr + 
    (os.result[[1]]$tfr.s * (os.result[[1]]$mab.s-mab.adj.var)^2 +
       os.result[[2]]$tfr.s * (os.result[[2]]$mab.s-mab.adj.var)^2 +
       os.result[[3]]$tfr.s * (os.result[[3]]$mab.s-mab.adj.var)^2 +
       os.result[[4]]$tfr.s * (os.result[[4]]$mab.s-mab.adj.var)^2) /
    i.tfr
  
  list(
    years = f.years,
    tfr.s=	i.tfr,
    tfr = tfr,
    mab.s = i.mab,
    mab = mab,
    r.bf = r.bf,
    vab.s = i.vab,
    vab = vab,
    tfr.adj.bf = (os.result[[1]]$tfr.adj.bf + 
                    os.result[[2]]$tfr.adj.bf +
                    os.result[[3]]$tfr.adj.bf + 
                    os.result[[4]]$tfr.adj.bf),  
    tfr.adj.var = tfr.adj.var,
    mab.adj.var = mab.adj.var,
    gamma = gamma,
    delta= delta,
    vab.adj.var = vab.adj.var
  )
}

###############################################################################
#### Input file: obj.APC 
#    [output of APC framework (R.ChengLin2010.function1.APC.r) or raw data]
#    Adjusted ASFR (named as "DEPEND") will be used
#### Parameters
#    mab: mean age at birth: fixed or missing (using mean mab of each year extimated in the function) 
#### Output: obj.KPS
#      obj.KPs$CTFR.all : all CTFR estimations
#      obj.KPs$CTFR.Raw : TFR
#      obj.KPs$CTFR.bf  : TFR with BF adjustment
#      obj.KPs$CTFR.bfs : CTFR.bf + kernel smoothing
#      obj.KPs$CTFR.KP  : TFR with var adjustment (KP)
#      obj.KPs$CTFR.KPS : CTFR.KP + kernel smoothing
#      obj.KPs$mab      : mab from input or estimation (mean mab of each year)
################
R.ChengLin2010.function2.KP.smooth<-function(obj.APC, mab = ""){
  
  ####
  if(sum(is.na(obj.APC$coef))>0){
    print("Non-complete data......return NA.")
    CTFR.all = matrix(NA,nrow=length(obj.APC$year1:obj.APC$year2),ncol=6)
    rownames(CTFR.all) <- paste(obj.APC$year1:obj.APC$year2)
    colnames(CTFR.all) <- c("YEAR","RAW","BF","BFS","KP","KPS")
    CTFR.all[,"YEAR"] <- obj.APC$year1:obj.APC$year2
    
    return(list(CTFR.all = CTFR.all,
                CTFR.Raw = NA, 
                CTFR.bf  = NA, 
                CTFR.bfs = NA, 
                CTFR.KP  = NA, 
                CTFR.KPS = NA, 
                mab= NA))
  }
  
  #print("KP smooth function 1")
  ## parameters
  age1 <- obj.APC$age1
  age2 <- obj.APC$age2
  year1 <- obj.APC$year1
  year2 <- obj.APC$year2
  
  ## get data
  use <- obj.APC$apcout
  ind <- use[,"YEAR"] >= year1 & use[,"YEAR"] <= year2
  use <- use[ind,]

  ## list of factors
  age <- sort(unique(use[,"AGE"]))
  period <- sort(unique(use[,"YEAR"]))
  cohort <- sort(unique(use[,"COHORT"]))
  
  
  temp <- matrix(0,nrow=length(age),ncol=length(period))
  rownames(temp) <- paste(age)
  colnames(temp) <- paste(period)
  
  for(i in 1:nrow(use)){
    temp[paste(use[i,"AGE"]),paste(use[i,"YEAR"])] <- use[i,"DEPEND"]
  }
  
  #### using KP S-Plus script
  TFR.s <- F.bf.calc(temp, age, period)
  
  ####
  if(mab <= 10)
    mab <- round(mean(TFR.s$mab,na.rm = T))
  #### kernel smoothing of BF
  TFR.s.bf.k <- ksmooth(TFR.s$years[-c(1,length(period))], TFR.s$tfr.adj.bf[-c(1,length(period))],
                        x.points=TFR.s$years[-c(1,length(period))],
                        kernel = "normal", bandwidth = 17)
  #### kernel smoothing of KP
  TFR.s.k <- ksmooth(TFR.s$years[-c(1,length(period))], TFR.s$tfr.adj.var[-c(1,length(period))],
                     x.points=TFR.s$years[-c(1,length(period))],
                     kernel = "normal", bandwidth = 17)
  #### 
  CTFR.Raw <- cbind(TFR.s$years - mab,TFR.s$tfr)
  colnames(CTFR.Raw) <- c("COHORT","CTFR")
  rownames(CTFR.Raw) <- paste(CTFR.Raw[,"COHORT"])
  
  CTFR.bf <- cbind(TFR.s$years - mab,TFR.s$tfr.adj.bf)
  colnames(CTFR.bf) <- c("COHORT","CTFR")
  rownames(CTFR.bf) <- paste(CTFR.bf[,"COHORT"])
  
  CTFR.bfs <- cbind(TFR.s.bf.k$x - mab,TFR.s.bf.k$y)
  colnames(CTFR.bfs) <- c("COHORT","CTFR")
  rownames(CTFR.bfs) <- paste(CTFR.bfs[,"COHORT"])
  
  ind <- !is.na(TFR.s$tfr.adj.var)
  CTFR.KP <- cbind(TFR.s$years[ind] - mab,TFR.s$tfr.adj.var[ind])
  colnames(CTFR.KP) <- c("COHORT","CTFR")
  rownames(CTFR.KP) <- paste(CTFR.KP[,"COHORT"])
  
  CTFR.KPS <- cbind(TFR.s.k$x - mab,TFR.s.k$y)
  colnames(CTFR.KPS) <- c("COHORT","CTFR")
  rownames(CTFR.KPS) <- paste(CTFR.KPS[,"COHORT"])
  
  CTFR.all <- cbind(TFR.s$years,
                    TFR.s$tfr,
                    TFR.s$tfr.adj.bf,
                    c(NA,TFR.s.bf.k$y,NA),
                    TFR.s$tfr.adj.var,
                    c(NA,TFR.s.k$y,NA))
  rownames(CTFR.all) <- paste(TFR.s$years)
  colnames(CTFR.all) <- c("YEAR","RAW","BF","BFS","KP","KPS")
  
  return(list(CTFR.all = CTFR.all,CTFR.Raw = CTFR.Raw, 
              CTFR.bf  = CTFR.bf, CTFR.bfs = CTFR.bfs, 
              CTFR.KP  = CTFR.KP, CTFR.KPS = CTFR.KPS, mab= TFR.s$mab))
  
  #plot(period,apply(temp,2,sum),type="l")
  #lines(TFR.s.k$x,TFR.s.k$y,col="red")
  #lines(TFR.s$years,TFR.s$tfr.adj.var,col="green")
}

