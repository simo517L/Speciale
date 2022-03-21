#liblocation = "/home/au591455/Rstuff/library"
liblocation = NULL
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
#setwd("/home/au591455/Rstuff/Results") 
setwd("C:/Users/simon/Desktop/TestR")
K_dist = function(M,k,p){
  dist_to_point = M[p,]
  dist_to_point = sort(dist_to_point)
  return(dist_to_point[k])
}

K_dist_neighborhood = function(M,k,p){
  K_dist_Result=K_dist(M,k,p)
  dist_to_point = M[p,]
  return(which(dist_to_point <= K_dist_Result))
}

lrd = function(M,k,p){
  PP_in_neighborhood = K_dist_neighborhood(M,k,p)
  Result= 0
  for (i in c(1:length(PP_in_neighborhood))){
    Result = Result + max(K_dist(M=M,k=k,p=PP_in_neighborhood[i]),M[p,PP_in_neighborhood[i]])
  }
  return(1/( Result/length(PP_in_neighborhood))) 
}


l_outlier_factor = function(M,k,p){
  PP_in_neighborhood = K_dist_neighborhood(M,k,p)
  Result= 0
  p_lrd=lrd(M,k,p) 
  for (i in c(1:length(PP_in_neighborhood))){
    Result = Result + lrd(M,k,PP_in_neighborhood[i]) / p_lrd
  }
  return( Result/length(PP_in_neighborhood)) 
}

outlier_factors_PP = function(X,k,nx,ny=ny,method=1,minpoints=20){
  M = distMppp(X,nx=nx,ny=ny,method=method,minpoints=minpoints)
  n = length(X)
  m = length(k)
  Result = matrix(0,nrow = n,ncol = m)
  colnames(Result) <- k
  for (i in  c(1:m)){
    for (j in c(1:n)){
      Result[j,i] = l_outlier_factor(M,k[i],j)
    }}
  return(Result)
}

r1 = rpoint(5,nsim=20)
r2 = rpoint(100)
Test_outlier_OF = function(Outlier,PPP,nx,ny=nx,method=1,minpoints=20,Kinterval){
  X = c(PPP,list(Outlier))
  n = length(PPP)
  Result = outlier_factors_PP(X=X,k=Kinterval,nx=nx,ny=ny,method=method,minpoints=minpoints)
  return(mean(Result[n+1] <= Result))
}

registerDoParallel(2)
nn =10000
mm=40
timeofpowertest  =c(1:6)
squares = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))


Data = rpoispp(100,nsim=20*nn)

Matern4a = rMaternI(100,r=0.05,nsim=nn)
for (i in c(1:mm)){
  
  grid1 = quadrats(Matern4a[[i]],nx=3,ny=4)
  splitMatern = split(Matern4a[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>3){
      count = count +1
    }
  }
  if (count < 3){
    Matern4a[[i]] = rMaternI(100,r=0.05)
    i = i-1
  }
}

Matern4b = rMaternI(100,r=0.02,nsim=nn)
for (i in c(1:mm)){
  grid1 = quadrats(Matern4b[[i]],nx=3,ny=4)
  splitMatern = split(Matern4b[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>3){
      count = count +1
    }
  }
  if (count < 3){
    Matern4b[[i]] = rMaternI(100,r=0.02)
    i = i-1
  }
}

Clust4a = rMatClust(100,scale=0.1,mu=1,nsim=nn)
for (i in c(1:mm)){
  grid1 = quadrats(Clust4a[[i]],nx=3,ny=4)
  splitMatern = split(Clust4a[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>3){
      count = count +1
    }
  }
  if (count < 3){
    Clust4a[[i]] = rMatClust(100,scale=0.1,mu=1)
    i = i-1
  }
}
Clust4b = rMatClust(100,scale=0.1,mu=4,nsim=nn)
for (i in c(1:mm)){
  grid1 = quadrats(Clust4b[[i]],nx=3,ny=4)
  splitMatern = split(Clust4b[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>3){
      count = count +1
    }
  }
  if (count < 3){
    Clust4b[[i]] = rMatClust(100,scale=0.1,mu=4)
    i = i-1
  }
}
Clust4c = rMatClust(100,scale=0.05,mu=4,nsim=nn)
for (i in c(1:mm)){
  grid1 = quadrats(Clust4c[[i]],nx=3,ny=4)
  splitMatern = split(Clust4c[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>3){
      count = count +1
    }
  }
  if (count < 3){
    Clust4c[[i]] = rMatClust(100,scale=0.05,mu=4)
    i = i-1
  }
}
poistest = rpoispp(100,nsim=nn)


logpath = "/home/au591455/Rstuff/Results/logOF.txt"
#logpath = "C:/Users/simon/Desktop/TestR/logOF.txt"
fileConn<-file(logpath)
writeLines("Start UP", fileConn)

powertest_OF_sq = function(Outlier,Data,name,n,m,squares,path_of_log){
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %:%
    foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
      QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 4)
      QQQ$p.value}
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))], nx=squares[[j]][1],ny=squares[[j]][2],minpoints =4)
        QQQ$p.value
      }
    ResultpowerM1  = cbind(ResultpowerM1,ResultpowerM1temp1)
    saveRDS(ResultpowerM1,file = name)
    tempC = readLines(path_of_log)
    procent = q*m/n
    writeLines(c(tempC,paste(name,":",procent," % done", Sys.time())),path_of_log)
  }
  
  
}

