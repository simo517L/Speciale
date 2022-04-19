
r = seq(0.05,0.3,0.05)
power_test = matrix(0,nrow = 10^4 , ncol = length(r))

for ( i in c(1:10^4)){
  Data = rpoispp(100,nsim=20)
  Outlier= rMaternI(100,r=0.05,nsim=1)
  for (j in  c(1:length(r))){
    QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=2,minpoints = 5,rinterval=c(0,r[j]))
    power_test[i,j] = QQQ$p.value < 0.05
  }
  if (i %% 100 == 0){
    print(i)
  }
}
colSums(power_test)/10^4

power_test2 = matrix(0,nrow = 10^4 , ncol = length(r))
squares = list(c(2,2),c(3,2),c(3,3),c(3,4),c(4,4))
for ( i in c(1:10^4)){
  Data = rpoispp(100,nsim=20)
  Outlier= rMatClust(100,scale=0.1,mu=1)
  for (j in  c(1:length(squares))){
    QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5,rinterval=c(0,r[j]))
    power_test2[i,j] = QQQ$p.value < 0.05
  }
  if (i %% 10 == 0){
    print(i)
  }
}

library(parallel)
library(foreach)
library(doParallel)
library(tictoc)



timeforloop = c(1:10)
for (i in c(1:10)){
  Data = rpoispp(100,nsim=20)
Outlier= rMatClust(100,scale=0.1,mu=1)
tic()
for (j in  c(1:length(squares))){
  QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
  QQQ$p.value < 0.05
}
T1 = toc()
timeforloop[i] = T1$toc-T1$tic
}

timeforreach = c(1:10)
registerDoParallel(10)
for (i in c(1:10)){
  Data = rpoispp(100,nsim=20)
  Outlier= rMatClust(100,scale=0.1,mu=1)
  tic()
foreach (j= c(1:length(squares)), .combine=c, .packages = c("spatstat")) %dopar%{
  QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
  QQQ$p.value < 0.05
}
T2 = toc()
timeforreach[i] = T2$toc-T2$tic
}

for (i in c(1:10)){
  Data = rpoispp(100,nsim=20)
  Outlier= rMatClust(100,scale=0.1,mu=1)
  tic()
  foreach (j= c(1:length(squares)), .combine=c, .packages = c("spatstat")) %dopar%{
    QQQ = OutlierPPP(Outlier = Outlier, PPP =  Data , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
  T2 = toc()
  timeforreach[i] = T2$toc-T2$tic
}

timeforloop2 = c(1:10)
for (q in c(1:10)){
  Data = rpoispp(100,nsim=20*5)
  Outlier= rMatClust(100,scale=0.1,mu=1,nsim=5)
  tic()
  
  for (i in c(1:5)){
    for  (j in c(1:length(squares))){
      QQQ = OutlierPPP(Outlier = Outlier[[i]], PPP = Data[c(i*20-19,i*20)] , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
      QQQ$p.value < 0.05
    }
    
  }
  T3 = toc()
  timeforloop2[q] = T3$toc-T3$tic
}


timeforsingle2 = c(1:10)
for (q in c(1:10)){
  Data = rpoispp(100,nsim=20*5)
  Outlier= rMatClust(100,scale=0.1,mu=1,nsim=5)
  tic()
for (i in c(1:5)){
  foreach (j= c(1:length(squares)), .combine=c, .packages = c("spatstat")) %dopar%{
    QQQ = OutlierPPP(Outlier = Outlier[[i]], PPP = Data[c(i*20-19,i*20)] , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }

}
  T3 = toc()
  timeforsingle2[q] = T3$toc-T3$tic
}

timefornest2 = c(1:10)
for (q in c(1:10)){
  Data = rpoispp(100,nsim=20*5)
  Outlier= rMatClust(100,scale=0.1,mu=1,nsim=5)
  tic()
    testM <- foreach (i= c(1:5), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        QQQ = OutlierPPP(Outlier = Outlier[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
        QQQ$p.value < 0.05
    
  }
  T3 = toc()
  timefornest2[q] = T3$toc-T3$tic
}


stopImplicitCluster()
mean(timeforloop2)
mean( timeforsingle2)
mean( timefornest2)



registerDoParallel(10)
n =1000
timeofpowertest  =c (1:5)
squares = list(c(2,2),c(3,2),c(3,3),c(3,4),c(4,4))

Data = rpoispp(100,nsim=20*n)
Outlier1= rMatClust(100,scale=0.1,mu=1,nsim=n)


Matern4a = rMaternI(100,r=0.05,nsim=n)
Matern4b = rMaternI(100,r=0.02,nsim=n)
Clust4a = rMatClust(100,scale=0.1,mu=1,nsim=n)
Clust4b = rMatClust(100,scale=0.1,mu=4,nsim=n)
Clust4c = rMatClust(100,scale=0.05,mu=4,nsim=n)

tic()
ResultpowerM1 <- foreach (i= c(1:n), .combine="cbind", .packages = c("spatstat")) %:%
  foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
    QQQ = OutlierPPP(Outlier = Matern4a[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
T1 = toc()
timeofpowertest[1] = T1$toc-T1$tic
print(1)
tic()
ResultpowerM2 <- foreach (i= c(1:n), .combine="cbind", .packages = c("spatstat")) %:%
  foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
    QQQ = OutlierPPP(Outlier = Matern4b[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
T1 = toc()
timeofpowertest[2] = T1$toc-T1$tic
print(2)
tic()
ResultpowerM3 <- foreach (i= c(1:n), .combine="cbind", .packages = c("spatstat")) %:%
  foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
    QQQ = OutlierPPP(Outlier = Clust4a[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
T1 = toc()
timeofpowertest[3] = T1$toc-T1$tic
print(3)
tic()
ResultpowerM4 <- foreach (i= c(1:n), .combine="cbind", .packages = c("spatstat")) %:%
  foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
    QQQ = OutlierPPP(Outlier = Clust4b[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
T1 = toc()
timeofpowertest[4] = T1$toc-T1$tic
print(4)
tic()
ResultpowerM5 <- foreach (i= c(1:n), .combine="cbind", .packages = c("spatstat")) %:%
  foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
    QQQ = OutlierPPP(Outlier = Clust4c[[i]], PPP = Data[c(i*20-19,i*20)], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5)
    QQQ$p.value < 0.05
  }
T1 = toc()
timeofpowertest[5] = T1$toc-T1$tic
stopImplicitCluster()

number_of_squares=c(4,6,9,12,16)
plot(number_of_squares,rowMeans(ResultpowerM2),type="l",lty=2,col=2,ylim = c(0,1),ylab = "Rate of rejection")
lines(number_of_squares,rowMeans(ResultpowerM3),lty=3,col=3)
lines(number_of_squares,rowMeans(ResultpowerM4),lty=4,col=4)
lines(number_of_squares,rowMeans(ResultpowerM5),lty=5,col=5)
legend(x=10,y=0.8 , legend = c("Matern2","Clustera","Clusterb","Clusterc"), lty=2:5,col = 2:5,cex = 0.55)



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

conf_Result2 = conf_of_power(ResultpowerM2,1000)
conf_Result3 = conf_of_power(ResultpowerM3,1000)
conf_Result4 = conf_of_power(ResultpowerM4,1000)
conf_Result5 = conf_of_power(ResultpowerM5,1000)
TestResult = list(ResultpowerM2,ResultpowerM3,ResultpowerM4,ResultpowerM5)
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

plot_of_power(Data=ResultpowerM2,Size = 1000,X= number_of_squares,conf = T,legendC = T)
plot_of_power(Data=TestResult,Size = 1000,X= number_of_squares,legendC = T)
legend("topright" , legend = c(1:n),col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))
plot(number_of_squares,rowMeans(ResultpowerM2),type="l",lty=2,col=2,ylim = c(0,1),ylab = "Rate of rejection")

lines(number_of_squares,rowMeans(ResultpowerM3),lty=3,col=3)
lines(number_of_squares,rowMeans(ResultpowerM4),lty=4,col=4)
lines(number_of_squares,rowMeans(ResultpowerM5),lty=5,col=5)


plot(L,p_value,type="l",lty=1,col=2,ylim = c(0,1),ylab = "Rate of rejection")
lines(L,upper_lim,lty=3)
lines(L,lower_lim,lty=3)


testconf = prop.test(rowSums(ResultpowerM2)[1],1000)
squares2  = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))


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

result_Mata_sumfunc_sq = readRDS(file = "PowerMaternA.Rdata")

number_of_squares = c(3,4,6,8,9,12) 
par(mfrow = c(3,2))
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

result_Mata_sumfunc_sq = result_Mata_sumfunc_sq[-6,]
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

result_Mata_sumfunc_sq = result_Mata_sumfunc_sq[-5,]
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

result_Mata_sumfunc_sq = result_Mata_sumfunc_sq[-4,]
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

result_Mata_sumfunc_sq = result_Mata_sumfunc_sq[-3,]
sum(is.nan(colMeans(result_Mata_sumfunc_sq)))

power_Mata_sunfunc_sq = result_Mata_sumfunc_sq[,!is.nan(colMeans(result_Mata_sumfunc_sq))] <=  0.05
sum(is.nan(colMeans(result_Mata_sumfunc_sq)) == T)
plot_of_power(Data=power_Mata_sunfunc_sq,title = "MaternA",Size = dim(power_Mata_sunfunc_sq)[2],X= number_of_squares[c(1,2)],conf = T,legendC = T)

result_Matb_sumfunc_sq = readRDS(file = "PowerMaternB.Rdata")
power_Matb_sunfunc_sq = result_Matb_sumfunc_sq[,is.nan(colMeans(result_Matb_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_Matb_sumfunc_sq)) == T)
plot_of_power(Data=power_Matb_sunfunc_sq,title = "MaternB",Size = dim(power_Matb_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

result_clustA_sumfunc_sq = readRDS(file = "PowerClusterA.Rdata")
power_clustA_sunfunc_sq = result_clustA_sumfunc_sq[,is.nan(colMeans(result_clustA_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustA_sumfunc_sq)) == T)
plot_of_power(Data=power_clustA_sunfunc_sq,title = "ClustA",Size = dim(power_clustA_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

result_clustB_sumfunc_sq = readRDS(file = "PowerClusterB.Rdata")
power_clustB_sunfunc_sq = result_clustB_sumfunc_sq[,is.nan(colMeans(result_clustB_sumfunc_sq)) != T] <=  0.05
sum(is.nan(colMeans(result_clustB_sumfunc_sq)) == T)
plot_of_power(Data=power_clustB_sunfunc_sq,title = "ClustB",Size = dim(power_clustB_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)


result_clustC_sumfunc_sq = readRDS(file = "PowerClusterC.Rdata")
power_clustC_sunfunc_sq = result_clustC_sumfunc_sq[,is.nan(colMeans(result_clustC_sumfunc_sq)) != T] <=  0.05

plot_of_power(Data=power_clustC_sunfunc_sq,title = "ClustC",Size = dim(power_clustC_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)


result_pois_sumfunc_sq = readRDS(file = "PowerpoisTest.Rdata")
power_pois_sunfunc_sq = result_pois_sumfunc_sq[,is.nan(colMeans(result_pois_sumfunc_sq)) != T] <=  0.05
plot_of_power(Data=power_pois_sunfunc_sq,title = "Pois",Size = dim(power_pois_sunfunc_sq)[2],X= number_of_squares,conf = T,legendC = T)

par(mfrow = c(1,1))

sum(is.nan(result_Mata_sumfunc_sq) == T)
sum(is.nan(result_Matb_sumfunc_sq) == T)

sum(is.nan(result_clustA_sumfunc_sq) == T)
sum(is.nan(result_clustB_sumfunc_sq) == T)
sum(is.nan(result_clustC_sumfunc_sq) == T)

sum(is.nan(result_pois_sumfunc_sq) == T)
