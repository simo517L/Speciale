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
conf_of_power=function(M,n){
  Result = rowSums(M)
  nn = length( Result)
  p_value = c(1:nn)
  upper_lim = c(1:nn)
  lower_lim = c(1:nn)
  for (i in c(1:nn)){
    TEMP = prop.test(Result[i],n)
    p_value[i] = TEMP$estimate
    upper_lim[i] = TEMP$conf.int[1]
    lower_lim[i] =  TEMP$conf.int[2]
  }
  return(data.frame(p_value=p_value,upper_lim=upper_lim,lower_lim=lower_lim))
}


plot_of_power = function(Data,Size,X,conf = FALSE,title="",legendC = FALSE,legendNames=NULL){
  if ( is.list(Data) & length(Data) > 1){
    n = length(Data)
    Results = list()
    for (i in c(1:n)){
      Results[[i]] = conf_of_power(Data[[i]],Size)
    }
    if (conf == FALSE){
      plot(X,Results[[1]]$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)
      for ( j in c(2:n)){
        lines(X,Results[[j]]$p_value,lty=j,col=j)
      }
    }     else{
      plot(X,Results[[1]]$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)
      lines(X,Results[[1]]$upper_lim,lty=2,col=1)
      lines(X,Results[[1]]$lower_lim,lty=2,col=1)
      for ( j in c(2:n)){
        lines(X,Results[[j]]$p_value,lty=1,col=j)
        lines(X,Results[[j]]$upper_lim,lty=2,col=j)
        lines(X,Results[[j]]$lower_lim,lty=2,col=j)
      }
      
    }} else if ( is.matrix(Data) & length(Data) > 1) {
      n = dim(Data)[2]
      Results = conf_of_power(Data,Size)
      
      if (conf == FALSE){
        plot(X,Results$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)
        
      }     else{
        plot(X,Results$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)
        lines(X,Results$upper_lim,lty=2,col=1)
        lines(X,Results$lower_lim,lty=2,col=1)
        
      }} else {
        n=1
        Results = conf_of_power(Data,Size)
        if (conf == FALSE){plot(X,Results$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)}
        else {plot(X,Results$p_value,type="l",lty=1,col=1,ylim = c(0,1),ylab = "Rate of rejection",main = title)
          lines(X,Results$upper_lim,lty=2,col=1)
          lines(X,Results$lower_lim,lty=2,col=1)
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



test_of_power = function(X, Outlier,Data,testrange=0.3,method){
  if (method == 1){
    r = seq(0.05,testrange,0.05)
    m = length(r)
    power_result  =c(1:m)
    for ( j in c(1:m)){
      QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=2,minpoints = 10,rinterval=c(0,r[j]))
      power_result[m] = QQQ$p.value < 0.05
    }
    
  } else if (method ==2){
    squares=X
    m = length(squares)
    power_result  =c(1:m)
    for ( j in c(1:m)){
      QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
      power_result[m] = QQQ$p.value < 0.05
    }
    
  }
  return(power_result)
}
setwd("C:/Users/simon/Desktop/TestR")
result_Mata_sumfunc_sq = readRDS(file = "PowerMaternA.Rdata")

number_of_squares = c(3,4,6,8,9,12) 
par(mfrow = c(3,2))
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

result_Mata_sumfunc_sq = result_Mata_sumfunc_sq[,!is.nan(colMeans(result_Mata_sumfunc_sq))]
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

power_Mata_sunfunc_sq = result_Mata_sumfunc_sq <=  0.05
sum(is.nan(colMeans(result_Mata_sumfunc_sq)) == T)
plot_of_power(Data=power_Mata_sunfunc_sq,title = "Matern: r=0.05",Size = dim(power_Mata_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

result_Matb_sumfunc_sq = readRDS(file = "PowerMaternB.Rdata")
power_Matb_sunfunc_sq = result_Matb_sumfunc_sq[,is.nan(colMeans(result_Matb_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_Matb_sumfunc_sq)) == T)
plot_of_power(Data=power_Matb_sunfunc_sq,title = "Matern: r=0.02",Size = dim(power_Matb_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

result_clustA_sumfunc_sq = readRDS(file = "PowerClusterA.Rdata")
power_clustA_sunfunc_sq = result_clustA_sumfunc_sq[,is.nan(colMeans(result_clustA_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustA_sumfunc_sq)) == T)
plot_of_power(Data=power_clustA_sunfunc_sq,title = "Clust: scale=0.1,mu=1",Size = dim(power_clustA_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

result_clustB_sumfunc_sq = readRDS(file = "PowerClusterB.Rdata")
power_clustB_sunfunc_sq = result_clustB_sumfunc_sq[,is.nan(colMeans(result_clustB_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustB_sumfunc_sq)) == T)
plot_of_power(Data=power_clustB_sunfunc_sq,title = "Clust: scale=0.1,mu=4",Size = dim(power_clustB_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)


result_clustC_sumfunc_sq = readRDS(file = "PowerClusterC.Rdata")
power_clustC_sunfunc_sq = result_clustC_sumfunc_sq[,is.nan(colMeans(result_clustC_sumfunc_sq)) != T] <=  0.05

plot_of_power(Data=power_clustC_sunfunc_sq,title = "Clust: scale=0.5,mu=4",Size = dim(power_clustC_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)


result_pois_sumfunc_sq = readRDS(file = "PowerpoisTest.Rdata")
power_pois_sunfunc_sq = result_pois_sumfunc_sq[,is.nan(colMeans(result_pois_sumfunc_sq)) != T] <=  0.05
plot_of_power(Data=power_pois_sunfunc_sq,title = "Poison",Size = dim(power_pois_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

par(mfrow = c(1,1))

