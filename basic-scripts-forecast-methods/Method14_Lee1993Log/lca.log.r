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

##  This file incorporates work covered by the following copyright and  
##  permission notice:  
##   
##      Copyright (c) 2003, Rob J Hyndman   
##      https://github.com/robjhyndman/demography   
##   
##      Permission to use, copy, modify, and/or distribute this software  
##      for any purpose with or without fee is hereby granted, provided  
##      that the above copyright notice and this permission notice appear  
##      in all copies.  
##   
##      THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL  
##      WARRANTIES WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED  
##      WARRANTIES OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE  
##      AUTHOR BE LIABLE FOR ANY SPECIAL, DIRECT, INDIRECT, OR  
##      CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS  
##      OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT,  
##      NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN  
##      CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.  
## 

########
## While using this script in 'validate-forecast-methods':
## In function lca(), we modified the format of input data.
## In function forecast.lca(), we modified the format of output data.
## In function findroot(), we change stop("Unable to find root") as return(NA) 
###############################################################################


## Last updated 29 May 2003 by RJH
## Added plot.lca and print.lca
## Changed way of limiting ages in lca
## Fixed bug which arose occasionally when using method "dt"
## Made "dt" the default as in original LC paper.

lca <-  function(data,expos,series="total",years=as.numeric(colnames(data)), ages=as.numeric(rownames(data)),
                 max.age=100, adjust=c("dt","dxt","e0","none"),
                 chooseperiod=FALSE, minperiod=20, breakmethod=c("bai","bms"),
                 scale=FALSE, restype=c("logrates","rates","deaths"), interpolate=FALSE,label="")
{

  adjust <- match.arg(adjust)
  restype <- match.arg(restype)
  breakmethod <- match.arg(breakmethod)

  startage <- min(ages)
  
  # Extract mortality rates and population numbers
  mx <- data
  pop <- expos
  
  # Truncate years
  startyear <- min(years)
  stopyear <- max(years)
  if(startyear > max(years) | stopyear < min(years))
    stop("Year not found")
  startyear <- max(startyear,min(years))
  if(!is.null(stopyear))
    stopyear <- min(stopyear,max(years))
  else
    stopyear <- max(years)
  id2 <- stats::na.omit(match(startyear:stopyear,years))
  
  mx <- mx[,id2]
  pop <- pop[,id2]
  year <- years[id2]
  deltat <- year[2]-year[1]
  ages <- ages
  n <- length(ages)
  m <- sum(id2>0)
  mx <- matrix(mx,nrow=n,ncol=m)
  #    n <- nrow(mx) # number of ages
  #    m <- ncol(mx) # number of years
  
  # Interpolate where rates are zero
  if(interpolate)
  {
    # Remove missing values
    mx[is.na(mx)] <- 0
    if(sum(abs(mx)<1e-9,na.rm=TRUE) > 0)
    {
      warning("Replacing zero values with estimates")
      for(i in 1:n)
        mx[i,] <- fill.zero(mx[i,])
    }
  }
  
  # Transpose data and get deaths and logrates
  mx <- t(mx)
  mx[mx==0] <- NA
  logrates <- log(mx)
  
  pop <- t(pop)
  deaths <- pop*mx
  
  # Do SVD
  ax <- apply(logrates,2,mean, na.rm=TRUE) # ax is mean of logrates by column
  if(sum(ax < -1e9) > 0)
    stop(sprintf("Some %s rates are zero.\n Try reducing the maximum age or setting interpolate=TRUE.", data$type))
  clogrates <- sweep(logrates,2,ax) # central log rates (with ax subtracted) (dimensions m*n)
  svd.mx <- svd(clogrates)
  
  # Extract first principal component
  sumv <- sum(svd.mx$v[,1])
  bx <- svd.mx$v[,1]/sumv
  kt <- svd.mx$d[1] * svd.mx$u[,1] * sumv
  
  # Adjust kt
  ktadj <- rep(0,m)
  logdeathsadj <- matrix(NA,n,m)
  z <- log(t(pop))+ax
  
  # Use regression to guess suitable range for root finding method
  x <- 1:m
  ktse <- stats::predict(stats::lm(kt ~ x),se.fit=TRUE)$se.fit
  ktse[is.na(ktse)] <- 1
  agegroup = ages[4]-ages[3]
  
  if (adjust=="dxt")
  {
    options(warn=-1) # Prevent warnings on non-integer population values
    for(i in 1:m)
    {
      y <- as.numeric(deaths[i,])
      zi <- as.numeric(z[,i])
      weight <- as.numeric(zi > -1e-8)  # Avoid -infinity due to zero population
      yearglm <- stats::glm(y ~ offset(zi)-1+bx, family=stats::poisson, weights=weight)
      ktadj[i] <- yearglm$coef[1]
      logdeathsadj[,i] <- z[,i]+bx*ktadj[i]
    }
    options(warn=0)
  }
  else if(adjust=="dt")
  {
    FUN <- function(p,Dt,bx,ax,popi) {Dt - sum(exp(p*bx+ax)*popi)}
    for(i in 1:m)
    {
      if(i==1)
        guess <- kt[1]
      else
        guess <- mean(c(ktadj[i-1],kt[i]))
      tryroot <- findroot(FUN, guess=guess, margin = 3*ktse[i],ax=ax,bx=bx,popi=pop[i,],Dt=sum(as.numeric(deaths[i,])))
      if(!is.na(tryroot)){
        ktadj[i] <- tryroot
      }
      if(is.na(tryroot)){
        ktadj[i] <- kt[i]
      }
      logdeathsadj[,i] <- z[,i]+bx*ktadj[i]
    }
  }
  else if(adjust=="e0")
  {
    e0 <- apply(mx,1,get.e0,agegroup=agegroup,sex=series,startage=startage)
    FUN2 <- function(p,e0i,ax,bx,agegroup,series,startage){e0i - estimate.e0(p,ax,bx,agegroup,series,startage)}
    for (i in 1:m)
    {
      if(i==1)
        guess <- kt[1]
      else
        guess <- mean(c(ktadj[i-1],kt[i]))
      ktadj[i] <- findroot(FUN2, guess=guess, margin = 3*ktse[i],e0i=e0[i],ax=ax,bx=bx,agegroup=agegroup,series=series,startage=startage)
    }
  }
  else if(adjust=="none")
    ktadj <- kt
  else
    stop("Unknown adjustment method")
  
  kt <- ktadj
  
  # Find linear section of kt and refit
  if(chooseperiod)
  {
    if(breakmethod=="bai")
    {
      x <- 1:m
      # Find breakpoints
      bp <- strucchange::breakpoints(kt ~ x)$breakpoints
      # Omit breakpoints less than minperiod from end
      bp <- bp[bp <= (m-minperiod)]
      bestbreak <- max(bp)
      return(lca(data,series,year[(bestbreak+1):m],ages=ages,max.age=max.age,
                 adjust=adjust,interpolate=interpolate,chooseperiod=FALSE,scale=scale))
    }
    else
    {
      RS <- devlin <- devadd <- numeric(m-2)
      for(i in 1:(m-2))
      {
        tmp <- lca(data,series,year[i:m],ages=ages,max.age=max.age,adjust=adjust,chooseperiod=FALSE,interpolate=interpolate,scale=scale)
        devlin[i] <- tmp$mdev[2]
        devadd[i] <- tmp$mdev[1]
        RS[i] <- (tmp$mdev[2]/tmp$mdev[1])
      }
      bestbreak <- order(RS[1:(m-minperiod)])[1]-1
      out <- lca(data,series,year[(bestbreak+1):m],ages=ages,max.age=max.age,
                 adjust=adjust,chooseperiod=FALSE,interpolate=interpolate,scale=scale)
      out$mdevs <- ts(cbind(devlin,devadd,RS),start=startyear,deltat=deltat)
      dimnames(out$mdevs)[[2]] <- c("Mean deviance total","Mean deviance base","Mean deviance ratio")
      return(out)
    }
  }
  
  # Estimate rates from fitted values and get residuals
  logfit <- fitmx(kt,ax,bx,transform=TRUE)
  if(restype=="logrates")
  {
    fit <- logfit
    res <- logrates - fit
  }
  else if(restype=="rates")
  {
    fit <- exp(logfit)
    res <- exp(logrates) - fit
  }
  else if(restype=="deaths")
  {
    fit <- exp(logfit)*pop
    res <- deaths - fit
  }
  residuals <- as.ts(ages,t(res),frequency=1/deltat,start=years[1],xname="Age",
                   yname=paste("Residuals", "fertility", "rate"))
  fitted <- as.ts(ages,t(fit),frequency=1/deltat,start=years[1],xname="Age",
                yname=paste("Fitted", "fertility", "rate"))
  
  names(ax) <- names(bx) <- ages
  
  # Rescaling bx, kt
  if(scale)
  {
    avdiffk <- -mean(diff(kt))
    bx <- bx*avdiffk
    kt <- kt/avdiffk
  }
  
  # Compute deviances
  deathsadjfit <- exp(logfit)*pop
  drift <- mean(diff(kt))
  ktlinfit <- mean(kt) + drift* (1:m - (m+1)/2)
  deathslinfit <- fitmx(ktlinfit,ax,bx,transform=FALSE)*pop
  dflogadd <- (m-2)*(n-1)
  mdevlogadd <- 2/dflogadd*sum(deaths*log(deaths/deathsadjfit)-(deaths-deathsadjfit))
  dfloglin <- (m-2)*n
  mdevloglin <- 2/dfloglin*sum(deaths*log(deaths/deathslinfit)-(deaths-deathslinfit))
  mdev <- c(mdevlogadd,mdevloglin)
  
  #Return
  output <- list(label=label,age=ages,year=year, mx=t(mx),
                 ax=ax, bx=bx, kt=ts(kt,start=startyear,deltat=deltat), residuals=residuals, fitted=fitted,
                 varprop=svd.mx$d[1]^2/sum(svd.mx$d^2), 
                 y=as.ts(ages,t(mx),start=years[1],frequency=1/deltat,xname="Age",
                       yname="Fertility"),
                 mdev=mdev)
  names(output)[4] <- series
  output$call <- match.call()
  names(output$mdev) <- c("Mean deviance base","Mean deviance total")
  output$adjust <- adjust
  output$type <- "fertility"
  return(structure(output,class="lca"))
}

bms <-  function(data,series=names(data$rate)[1],years=data$year, ages=data$age,
                 max.age=100, minperiod=20, breakmethod=c("bms","bai"), scale=FALSE, restype=c("logrates","rates","deaths"),
                 interpolate=FALSE)
{
  restype <- match.arg(restype)
  breakmethod <- match.arg(breakmethod)
  out <- lca(data,series=series, years=years, ages=ages, max.age=max.age, adjust="dxt",
             chooseperiod=TRUE, minperiod=minperiod, scale=scale, restype=restype, breakmethod=breakmethod,
             interpolate=interpolate)
  out$call <- match.call()
  return(out)
}

estimate.e0 <- function(kt,ax,bx,agegroup,series,startage=0)
{
  if(length(kt)>1)
    stop("Length of kt greater than 1")
  mx <- c(fitmx(kt,ax,bx))
  return(get.e0(mx,agegroup,series,startage=startage))
}

fitmx <- function (kt,ax,bx,transform=FALSE)
{
  # Derives mortality rates from kt mortality index,
  # per Lee-Carter method
  clogratesfit <- outer(kt, bx)
  logratesfit <- sweep(clogratesfit,2,ax,"+")
  if(transform)
    return(logratesfit)
  else
    return(exp(logratesfit))
}

plot.lca <- function(x,...)
{
  x$basis <- cbind(x$ax,x$bx)
  x$coeff <- cbind(rep(1,length(x$kt)),x$kt)
  colnames(x$basis) <- c("mean","bx")
  if(x$adjust != "none")
    xlab <- "kt (adjusted)"
  else
    xlab <- "kt"
  ftsa::plot.ftsm(x,1,"Age","bx","Year",xlab,mean.lab="ax",...)
}


print.lca <- function(x,...)
{
  cat("Lee-Carter analysis\n")
  cat(paste("\nCall:",deparse(x$call),"\n"))
  cat(paste("\nAdjustment method:",x$adjust))
  cat(paste("\nRegion:"),x$label)
  cat(paste("\nYears in fit:",min(x$year),"-",max(x$year)))
  cat(paste("\nAges in fit:",min(x$age),"-",max(x$age),"\n"))
  cat(paste("\nPercentage variation explained: ",round(x$varprop*100,1),"%\n",sep=""))
}

summary.lca <- function(object,...)
{
  print(object)
  
  cat(sprintf("\nERROR MEASURES BASED ON %s RATES\n", toupper(object$type)))
  printout(fdmMISE(object[[4]],exp(object$fitted$y),age=object$y$x,years=object$year))
  
  cat(sprintf("\nERROR MEASURES BASED ON LOG %s RATES\n", toupper(object$type)))
  printout(fdmMISE(log(object[[4]]),object$fitted$y,age=object$y$x,years=object$year))
}

printout <- function(output)
{
  junk1 <- cbind(output$ME,output$MSE,output$MPE,output$MAPE)
  rownames(junk1) <- output$age
  colnames(junk1) <- c("ME","MSE","MPE","MAPE")
  junk2 <- cbind(output$MIE,output$MISE,output$MIPE,output$MIAPE)
  rownames(junk2) = output$year
  colnames(junk2) = c("IE","ISE","IPE","IAPE")
  cat(paste("\nAverages across ages:\n"))
  print(round(apply(junk1,2,mean),5))
  cat(paste("\nAverages across years:\n"))
  print(round(apply(junk2,2,mean),5))
  cat("\n")
}

# Function performs predictions of k and life expectancy based on leecarter results (in lcaout)

forecast.lca <- function(object, h=60, se=c("innovdrift","innovonly"),method=c("arma","rw"), jumpchoice=c("fit","actual"), level=80, ...)
{
  se <- match.arg(se)
  jumpchoice <- match.arg(jumpchoice)
  
  # Section 1 Read in data from object
  jumpyear <- max(object$year)
  nyears <- length(object$year)
  nages <- length(object$age)
  
  # Find jumprates
  if(jumpchoice=="actual")
    jumprates <- object[[4]][,nyears]
  else if(jumpchoice=="fit")
    jumprates <- exp(object$ax + object$bx*object$kt[nyears])
  else
    stop(paste("Unknown jump choice:",jumpchoice))
  object$kt <- object$kt - object$kt[nyears]
  
  if(method=="rw"){
    # Time series estimation of kt as Random walk with drift
    fit <- forecast::rwf(object$kt, drift=TRUE)
    kt.drift <- fit$model$drift
    sec <- fit$model$drift.se
    see <- fit$model$sd
    
    # Project kt
    x <- 1:h
    zval <- stats::qnorm(0.5 + 0.005*level)
    kt.forecast <- object$kt[nyears] + (x * kt.drift)
    
    # Calculate standard errors of forecast kt
    if (se=="innovdrift")
      kt.stderr <- sqrt(x*(see^2) + (x*sec)^2)
    else if(se=="innovonly")
      kt.stderr <- sqrt(x*(see^2))
    kt.lo.forecast <- kt.forecast - (zval*kt.stderr)
    kt.hi.forecast <- kt.forecast + (zval*kt.stderr)
    kt.f <- data.frame(kt.forecast,kt.lo.forecast,kt.hi.forecast)
  }
  if(method=="arma"){
    obj <- auto.arima(object$kt,max.p = 5,max.d = 2,max.q=2)
    kt.f <- as.data.frame(forecast(obj,h = h,level = level))
    #print("arma")
  }
  
  names(kt.f) <- c("kt forecast","kt lower","kt upper")
  deltat <- object$year[2] - object$year[1]
  kt.f <- ts(kt.f,start=object$year[nyears]+deltat,deltat=deltat)
  
  # Calculate expected life and mx forecasts
  mx.forecast <- mx.lo.forecast <- mx.hi.forecast <- matrix(0,nrow=nages,ncol=h)
  logjumprates <- log(jumprates)
  for (cnt in 1:h)
  {
    mx.forecast[,cnt] <- fitmx(kt.f[cnt,1], logjumprates, object$bx)
    mx.lo.forecast[,cnt] <- fitmx(kt.f[cnt,2], logjumprates, object$bx)
    mx.hi.forecast[,cnt] <- fitmx(kt.f[cnt,3], logjumprates, object$bx)
  }

  output = list(label=object$label,age=object$age,year=object$year[nyears] + c(1:h)*deltat,
                rate=list(forecast=mx.forecast,lower=mx.lo.forecast,upper=mx.hi.forecast),
                fitted=object$fitted,
                kt.f=structure(list(mean=kt.f[,1],lower=kt.f[,2],upper=kt.f[,3],level=level,x=object$kt,
                                    method="Random walk with drift"),class="forecast"),
                type = object$type,lambda=0)
  names(output$rate)[1] = names(object)[4]
  output$model <- object
  output$model$jumpchoice <- jumpchoice
  output$model$jumprates <- jumprates
  output$call <- match.call()
  output$name <- names(object)[4]
  return(structure(output,class=c("fmforecast","demogdata")))
}

fitted.lca <- function(object,...)
{
  object$fitted
}

residuals.lca <- function(object,...)
{
  return(structure(list(x=object$year,y=object$age,z=t(object$residuals$y)),class="fmres"))
}

findroot <- function(FUN,guess,margin,try=1,...)
{
  # First try in successively larger intervals around best guess
  for(i in 1:5)
  {
    rooti <- try(stats::uniroot(FUN,interval=guess+i*margin/3*c(-1,1),...),silent=TRUE)
    if(class(rooti) != "try-error")
      return(rooti$root)
  }
  # No luck. Try really big intervals
  rooti <- try(stats::uniroot(FUN,interval=guess+10*margin*c(-1,1),...),silent=TRUE)
  if(class(rooti) != "try-error")
    return(rooti$root)
  
  # Still no luck. Try guessing root using quadratic approximation
  if(try<3)
  {
    root <- try(quadroot(FUN,guess,10*margin,...),silent=TRUE)
    if(class(root)!="try-error")
      return(findroot(FUN,root,margin,try+1,...))
    root <- try(quadroot(FUN,guess,20*margin,...),silent=TRUE)
    if(class(root)!="try-error")
      return(findroot(FUN,root,margin,try+1,...))
  }
  
  # Finally try optimization
  root <- try(newroot(FUN,guess,...),silent=TRUE)
  if(class(root)!="try-error")
    return(root)
  else
    root <- try(newroot(FUN,guess-margin,...),silent=TRUE)
  if(class(root)!="try-error")
    return(root)
  else
    root <- try(newroot(FUN,guess+margin,...),silent=TRUE)
  if(class(root)!="try-error")
    return(root)
  else
    return(NA)
}

quadroot <- function(FUN,guess,margin,...)
{
  x1 <- guess-margin
  x2 <- guess+margin
  y1 <- FUN(x1,...)
  y2 <- FUN(x2,...)
  y0 <- FUN(guess,...)
  if(is.na(y1) | is.na(y2) | is.na(y0))
    stop("Function not defined on interval")
  b <- 0.5*(y2-y1)/margin
  a <- (0.5*(y1+y2)-y0)/(margin^2)
  tmp <- b^2 - 4*a*y0
  if(tmp < 0)
    stop("No real root")
  tmp <- sqrt(tmp)
  r1 <- 0.5*(tmp-b)/a
  r2 <- 0.5*(-tmp-b)/a
  if(abs(r1) < abs(r2))
    root <- guess+r1
  else
    root <- guess+r2
  return(root)
}

# Try finding root using minimization
newroot <- function(FUN,guess,...)
{
  fred <- function(x,...){(FUN(x,...)^2)}
  junk <- stats::nlm(fred,guess,...)
  if(abs(junk$minimum)/fred(guess,...) > 1e-6)
    warning("No root exists. Returning closest")
  return(junk$estimate)
}


# Replace zeros with interpolated values
fill.zero <- function(x,method="constant")
{
  tt <- 1:length(x)
  zeros <- abs(x) < 1e-9
  xx <- x[!zeros]
  tt <- tt[!zeros]
  x <- stats::approx(tt,xx,1:length(x),method=method,f=0.5,rule=2)
  return(x$y)
}
