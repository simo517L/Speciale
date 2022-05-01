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
nn =10000
squares = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))
h=0.05
lambda = 100
-log(1-lambda*pi*h^2)/(pi*h^2)
1-exp(-106*pi)

Data = rpoispp(100,nsim=2*nn)
saveRDS(Data,file = "DataPPP.Rdata")

Matern4a = rMaternII(196,r=0.05,nsim=nn)
for (i in c(1:nn)){
  
  grid1 = quadrats(Matern4a[[i]],nx=3,ny=4)
  splitMatern = split(Matern4a[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>4){
      count = count +1
    }
  }
  if (count < 2){
    Matern4a[[i]] = rMaternI(100,r=0.05)
    i = i-1
  }
}
saveRDS(Matern4a,file = "Matern_a.Rdata")
Matern4b = rMaternII(106,r=0.02,nsim=nn)
for (i in c(1:nn)){
  grid1 = quadrats(Matern4b[[i]],nx=3,ny=4)
  splitMatern = split(Matern4b[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>4){
      count = count +1
    }
  }
  if (count < 3){
    Matern4b[[i]] = rMaternI(107,r=0.02)
    i = i-1
  }
}
saveRDS(Matern4b,file = "Matern_b.Rdata")
Clust4a = rMatClust(100,scale=0.1,mu=1,nsim=nn)
for (i in c(1:nn)){
  grid1 = quadrats(Clust4a[[i]],nx=3,ny=4)
  splitMatern = split(Clust4a[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>4){
      count = count +1
    }
  }
  if (count < 2){
    Clust4a[[i]] = rMatClust(100,scale=0.1,mu=1)
    i = i-1
  }
}
saveRDS(Clust4a,file = "Clust_a.Rdata")
Clust4b = rMatClust(25,scale=0.1,mu=4,nsim=nn)
for (i in c(1:nn)){
  grid1 = quadrats(Clust4b[[i]],nx=3,ny=4)
  splitMatern = split(Clust4b[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>4){
      count = count +1
    }
  }
  if (count < 3){
    Clust4b[[i]] = rMatClust(25,scale=0.1,mu=4)
    i = i-1
  }
}
saveRDS(Clust4b,file = "Clust_b.Rdata")
Clust4c = rMatClust(25,scale=0.05,mu=4,nsim=nn)
for (i in c(1:mm)){
  grid1 = quadrats(Clust4c[[i]],nx=3,ny=4)
  splitMatern = split(Clust4c[[i]],f=grid1)
  count = 0
  for (j in c(1:length(splitMatern))){
    if (splitMatern[[j]]$n>4){
      count = count +1
    }
  }
  if (count < 3){
    Clust4c[[i]] = rMatClust(25,scale=0.05,mu=4)
    i = i-1
  }
}
saveRDS(Clust4c,file = "Clust_c.Rdata")
poistest = rpoispp(100,nsim=nn)
saveRDS(poistest,file = "poisPPP.Rdata")


