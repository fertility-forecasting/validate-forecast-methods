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



##################################
APE_ks_test_grid_plot <- function(input,method.name="",ref.col="",legend.txt=NA,legend.col=NA, clust=F, ord=T, crossover=T){
  
  if(length(ref.col)<ncol(input)){
    ref.col <- rep(1,ncol(input))
  }
  zones=matrix(c(1,2), ncol=1, byrow=TRUE)
  layout(zones, widths=c(1), heights=c(8.7/10,1.3/10))
  
  if(clust == T){
    dis <- dist(input,method="man")
    fit <- hclust(dis, method = "single")
    #plot(fit)
    input <- input[fit$order,fit$order]
    method.name <- method.name[fit$order]
    ref.col <- ref.col[fit$order]
  }
  
  if(ord == T){
    num <- apply(input,1,function(x){sum(x<0.05)})
    for(i in 1:nrow(input)){
      num[i] <- sum(input[i,]<0.05 & input[,i]>=0.05)
    }
    reord <- order(num - 0.0001*apply(input,1,function(x){sum(x)}))
    input <- input[reord,reord]
    method.name <- method.name[reord]
    ref.col <- ref.col[reord]
  }
  input[!upper.tri(input) & !lower.tri(input)] <- 1.1
  myColor <- c("navy","grey60","grey99")
  myBreaks <- c(0, 0.05,1, 1.1)
  
  
  par(mar=c(0.1,15,13,0.9))
  
  image(1:ncol(input), 1:nrow(input), t(input)[nrow(input):1,], col=myColor, breaks=myBreaks, axes = FALSE, xlab="", ylab="")
  
  mtext(side = 3, at = 1:nrow(input), text = rownames(input)[nrow(input):1], col = ref.col[nrow(input):1], cex = 1.2, line = 1.2,las=3)
  mtext(side = 2, at = 1:nrow(input), text = rownames(input), col = ref.col, cex = 1.2, line = 1, par(las=1))
  if(all(!is.na(method.name))){
    mtext(side = 2, at = 1:nrow(input), text = method.name, col = ref.col, cex = 1.2, line = 1, par(las=1))
  }
  
  if(all(!is.na(legend.txt))){
    legend(-11.5,29.4, legend=legend.txt,text.col = legend.col,xjust = 0,xpd=T,text.width = 8, x.intersp=-0.2,cex = 1)
  }
  
  grid(nrow(input),ncol(input),col="grey95",lty=1)
  
  if(crossover==T){
    for(i in 1:nrow(input)){
      for(j in 1:ncol(input)){
        if(input[i,j] < 0.05 & input[j,i] < 0.05){
          text(nrow(input)+1-i,j,"X",col=myColor[3])
          text(nrow(input)+1-j,i,"X",col=myColor[3])
        }
      }
    }
  }
  ####
  scl <- nrow(input)/20
  ####
  text(-0.0,23.72*scl+0.5,"Method B:",xpd=T,cex=1.35,srt=90,font=2)
  text(-2.5,21*scl+0.5,"Method A:",xpd=T,cex=1.35,srt=0,font=2)
  text(-11*scl,11*scl,"Improving performance",xpd=T,srt=90,cex=1.5,font=2.5)
  arrows(-10*scl,19*scl,-10*scl,2*scl,xpd=T,lwd = 3,angle = 25,code=1)
  ####
  par(mar=c(0,15,0,0))
  plot.new()
  plot.window( xlim=c(1,ncol(input)+2), ylim=c(0,2) )
  text(x=c(1,9.8,16.4)*scl+0.5,c(0.74),c("Non significant","Significant","Crossover"),cex=1.2,pos=4)
  rect(0.14*scl+0.3,  0.5, 1.21*scl+0.3, 1.0,col=myColor[2],border = "white",lwd=0.1)
  #rect(7.2, 0.5, 8.2,1.0,col=myColor[2],border = "white",lwd=0.1)
  text(16.2*scl+0.5, 0.75,"X",col="black",cex=1.2,font=1.2)
  rect(8.92*scl+0.5, 0.5, 9.99*scl+0.5,1.0,col=myColor[1],border = "white",lwd=0.1)
  
  text(x=11.13*scl+0.5,y = 1.4,"Method A (row) performs strictly better than method B (column):")
  ####
  
  return(input)
}

##############################
##############################
ref.col <- rep("black" ,20)
##########