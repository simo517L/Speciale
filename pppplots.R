
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

setwd("C:/Users/simon/Desktop/TestR") 

nn = 3


MaternshowA= rMaternI(100,r=0.05,nsim=nn)

MaternshowB = rMaternI(100,r=0.02,nsim=nn)

ClustshowA = rMatClust(100,scale=0.1,mu=1,nsim=nn)

ClustshowB = rMatClust(100,scale=0.1,mu=4,nsim=nn)

ClustshowC = rMatClust(100,scale=0.05,mu=4,nsim=nn)

poispppShow = rpoispp(100,nsim=nn)

LLL = c(MaternshowA,MaternshowB,ClustshowA,ClustshowB,ClustshowC,poispppShow)

LLL
par(mfrow=c(3,6),mai = c(0.1, 0.1, 0.1, 0.1))
NAMESC = sample(paste("exp", c(1:18)))

for (i in c(1:length(LLL))) {
  temp = LLL[[i]]
  plot(temp, main = NAMESC[i])
}

MaternA = c(12,6,9)

MaternB = c(8,18,16)

ClustA = c(7,15,17)

ClustB = c(14,3,2)

ClustC = c(10,13,1)

pois =  c(5,11,4)

set.seed(1)
LLL2 = sample(LLL)
set.seed(1)
NAMESC2 = sample(NAMESC)
par(mfrow=c(3,3),mai = c(0.1, 0.1, 0.1, 0.1))
n = length(LLL)
for (i in c(1:9)) {
  temp = LLL2[[i]]
  plot(temp, main = NAMESC2[i])
}

for (i in c(10:18)) {
  temp = LLL2[[i]]
  plot(temp, main = NAMESC2[i])
}













load("megacaryocytes.rda")

load("intermembrane_particles.rda")
plot(megacaryocytes[1:10])
plot(megacaryocytes[11:20])
intermembrane_particles
plot(intermembrane_particles[1:10])
plot(intermembrane_particles[34:44])
plot(intermembrane_particles[48:58])
