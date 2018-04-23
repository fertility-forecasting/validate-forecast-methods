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

#################################################################################
#### load data and compute APE
#################################################################################
########
CV_APE_threshold_plot <- function(input,method.name=NA,plot.file="test.pdf",main="",xlim=c(-43,85),step=c(10,20),col=""){
  
  ####
  scl <- length(input)/20
  print(scl)
  ####
  if(any(col=="")){
    col=rep(1,length(input))
  }
  
  ####
  pdf(file=plot.file,height=8,width=8,family="Helvetica",pointsize=14)
  
  par(mfrow=c(1,1),las=1,mar=c(1.1,0.5,0.6,0.1))
  
  plot(x=-100,y=-100,xlim=xlim,ylim=c(0.5,15.5),xlab="",ylab="",main=main,axes=FALSE,cex.main=1.2)
  axis(side=1,at=seq(0,xlim[2],step[1]),labels=FALSE,lwd=1,pos=0.5,cex.axis=1.2)
  axis(side=1,at=seq(0,xlim[2],step[2]),labels=TRUE, lwd=3,pos=0.5,cex.axis=1.2)
  
  text(x=min(xlim)*1.15,y=15.3,"% threshold APE:",pos=4,cex=1.2,font=2)
  legend(x=xlim[1]/20,y=15.95,c( "50%   ","80%   ","90%   ","95%   "),horiz=TRUE,lwd=2,x.intersp =0.5,
         col=c("black","black","black","black"),lty=c(NA,NA,NA,NA),
         pch=c(18,17,16,15),cex=1.2,bty="n")
  
  text(max(xlim)*0.9,0.8,"APE in %",cex=1.2)
  yUp <- 0.68/scl
  
  ######################  

  if(any(is.na(method.name))){
    method.name <- names(input)
  }
  
  
  for(n in 1:length(input)){
    text(x=min(xlim)*1.15,y=0.5+(yUp*n),method.name[n],font=2,col=col[3],pos=4,cex=1.2)
    segments(y0=0.5+(yUp*n),y1=0.5+(yUp*n),x0=quantile(input[[n]]*100,prob=c(0),na.rm=TRUE),x1=quantile(input[[n]]*100,prob=c(0.95),na.rm=TRUE))
    segments(y0=0.4+(yUp*n),y1=0.6+(yUp*n),x0=quantile(input[[n]]*100,prob=c(0),na.rm=TRUE),x1=quantile(input[[n]]*100,prob=c(0),na.rm=TRUE))
    points(y=0.5+(yUp*n),x=quantile(input[[n]]*100,prob=c(0.95),na.rm=TRUE),pch=15,col=col[3],cex=1.5)
    points(y=0.5+(yUp*n),x=quantile(input[[n]]*100,prob=c(0.90),na.rm=TRUE),pch=16,col=col[3],cex=1.5)
    points(y=0.5+(yUp*n),x=quantile(input[[n]]*100,prob=c(0.80),na.rm=TRUE),pch=17,col=col[3],cex=1.5)
    points(y=0.5+(yUp*n),x=quantile(input[[n]]*100,prob=c(0.50),na.rm=TRUE),pch=18,col=col[3],cex=1.5)
  }

 
  dev.off()
  
}
##################################
