liblocation=NULL
library("spatstat.data",lib.loc=liblocation )
library("spatstat.geom",lib.loc=liblocation )
library("spatstat.random",lib.loc=liblocation )
library("spatstat.core",lib.loc=liblocation )
library("spatstat.linnet",lib.loc=liblocation )
library("spatstat",lib.loc=liblocation )
library(parallel)
library("codetools",lib.loc=liblocation )
library("iterators",lib.loc=liblocation )
library("tictoc",lib.loc=liblocation )
library("foreach",lib.loc=liblocation )
library("doParallel",lib.loc=liblocation )

# This function used prop.test to make a confidens interval for our data
conf_of_power=function(M,n){
  Result = rowSums(M)
  nn = length( Result)
  p_value = c(1:nn)
  upper_lim = c(1:nn)
  lower_lim = c(1:nn)
  for (i in c(1:nn)){
    TEMP = prop.test(Result[i], n, correct = FALSE)
    p_value[i] = TEMP$estimate
    upper_lim[i] = TEMP$conf.int[1]
    lower_lim[i] =  TEMP$conf.int[2]
  }
  return(data.frame(p_value=p_value,upper_lim=upper_lim,lower_lim=lower_lim))
}

# This function plots the data
plot_of_power = function(Data,Size,X,conf = FALSE,title="",legendC = FALSE,legendNames=NULL,yval=NULL){
  if (!is.null(yval)){
    a=yval[1]
    b=yval[2]
  }
  if ( is.list(Data) & length(Data) > 1){
    n = length(Data)
    Rtemp=c(1:n)
    Rlow=c(1:n)
    Rup=c(1:n)
    for (i in c(1:n)){
      Results = conf_of_power(Data[[i]],Size)
      Rtemp[i] = Results$p_value
      Rlow[i] = Results$lower_lim
      Rup[i] = Results$upper_lim
    }
    if (conf == FALSE){
      
      if (is.null(yval)){
        a = min(Rtemp) -0.05
        b = max(Rtemp) +0.05
      }
      
      plot(X,Rtemp[1],type="p",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)
      for ( j in c(2:n)){
        points(X,Rtemp[j],lty=j,col=j)
      }
    }     else{
      if (is.null(yval)){
        a = min(Rlow) -0.05
        b = max(Rup) +0.05
      }
      plot(X,Rtemp[1],type="p",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)
      arrows(X, Rlow[1],X, Rup[1], length=0.01, angle=90, code=3) 
      #points(X,Rup[1],lty=2,col=1)
      #points(X,Rlow[1],lty=2,col=1)
      for ( j in c(2:n)){
        points(X,Rtemp[j],lty=1,col=j)
        arrows(X, Rlow[j],X, Rup[j], length=0.01, angle=90, code=3)
        #points(X,Rup[j],lty=2,col=j)
        #points(X,Rlow[j],lty=2,col=j)
      }
      
    }} else if ( is.matrix(Data) & length(Data) > 1) {
      n = dim(Data)[2]
      Results = conf_of_power(Data,Size)
      
      if (conf == FALSE){
        if (is.null(yval)){
          a = min(Results$p_value) -0.1
          b = max(Results$p_value) +0.1
        }
        plot(X,Results$p_value,type="p",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)
        
      }     else{
        if (is.null(yval)){
          a = min(Results$lower_lim) -0.1
          b = max(Results$upper_lim) +0.1
        }
        plot(X,Results$p_value,type="p",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)
        arrows(X, Results$lower_lim,X, Results$upper_lim, length=0.01, angle=90, code=3)
        #points(X,Results$upper_lim,lty=2,col=1)
        #points(X,Results$lower_lim,lty=2,col=1)
        
      }} else {
        n=1
        Results = conf_of_power(Data,Size)
        if (conf == FALSE){
          if (is.null(yval)){
            a = min(Results$p_value) -0.1
            b = max(Results$p_value) +0.1
          }
          plot(X,Results$p_value,type="p",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)}
        else {
          if (is.null(yval)){
            a = min(Results$lower_lim) -0.1
            b = max(Results$upper_lim) +0.1
          }
          plot(X,Results$p_value,type="l",lty=1,col=1,ylim = c(a,b),ylab = "Rate of rejection",main = title)
          arrows(X, Results$lower_lim,X, Results$upper_lim, length=0.01, angle=90, code=3)
          #points(X,Results$upper_lim,lty=2,col=1)
          #points(X,Results$lower_lim,lty=2,col=1)
        }
      }
  if (legendC & n>1){
    if (is.null(legendNames)){#https://www.statology.org/legend-outside-plot-r/
      if (is.null(names(Data))){
        legend("topright" , legend = c(1:n),col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))
      } else{
        legend("topright" , legend = names(Data),col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))
      }
    } else{legend("topright" , legend = legendNames,col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))}
  }
}




path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path)  
# then we load the the results for the permutation test and plot then using the plot_of_power  function 
number_of_squares = c(4,6,9,12) 
par(mfrow = c(3,2))
result_Mata_sumfunc_sq = readRDS(file = "PowerMaternA.Rdata")

power_Mata_sunfunc_sq = result_Mata_sumfunc_sq <=  0.05
sum(is.nan(colMeans(result_Mata_sumfunc_sq)) == T)
plot_of_power(Data=power_Mata_sunfunc_sq,yval=c(0,1),title = "Matern: r=0.05",Size = dim(power_Mata_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)

result_Matb_sumfunc_sq = readRDS(file = "PowerMaternB.Rdata")
power_Matb_sunfunc_sq = result_Matb_sumfunc_sq[,is.nan(colMeans(result_Matb_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_Matb_sumfunc_sq)) == T)
plot_of_power(Data=power_Matb_sunfunc_sq,yval=c(0,1),title = "Matern: r=0.02",Size = dim(power_Matb_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)

result_clustA_sumfunc_sq = readRDS(file = "PowerClusterA.Rdata")
power_clustA_sunfunc_sq = result_clustA_sumfunc_sq[,is.nan(colMeans(result_clustA_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustA_sumfunc_sq)) == T)
plot_of_power(Data=power_clustA_sunfunc_sq,yval=c(0,1),title = "Clust: scale=0.1,mu=1",Size = dim(power_clustA_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)

result_clustB_sumfunc_sq = readRDS(file = "PowerClusterB.Rdata")
power_clustB_sunfunc_sq = result_clustB_sumfunc_sq[,is.nan(colMeans(result_clustB_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustB_sumfunc_sq)) == T)
plot_of_power(Data=power_clustB_sunfunc_sq,yval=c(0,1),title = "Clust: scale=0.1,mu=4",Size = dim(power_clustB_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)


result_clustC_sumfunc_sq = readRDS(file = "PowerClusterC.Rdata")
power_clustC_sunfunc_sq = result_clustC_sumfunc_sq[,is.nan(colMeans(result_clustC_sumfunc_sq)) != T] <=  0.05

plot_of_power(Data=power_clustC_sunfunc_sq,yval=c(0,1),title = "Clust: scale=0.5,mu=4",Size = dim(power_clustC_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)


result_pois_sumfunc_sq = readRDS(file = "PowerpoisTest.Rdata")
power_pois_sunfunc_sq = result_pois_sumfunc_sq[,is.nan(colMeans(result_pois_sumfunc_sq)) != T] <=  0.05
plot_of_power(Data=power_pois_sunfunc_sq,yval=c(0,1),title = "Poison",Size = dim(power_pois_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = F)

par(mfrow = c(1,1))

