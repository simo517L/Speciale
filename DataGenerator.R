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
path = path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember yo set your own path 
setwd(path) 

# In this file we will be generation data
# We start by generation the point patterns we will need for our simulations studies.
nn =10000
squares = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))



Data = rpoispp(100,nsim=2*nn)
saveRDS(Data,file = "DataPPP.Rdata")
# We use a loop to ensure each pattern has at least 2 squares with more then 4 point in them.
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

#Making point pattern from membrane data using thinning

load(file = "intermembrane_particles.rda") # We did not include this intermembrane_particles.rda with the code. 
tempv = c(1:length(intermembrane_particles$pattern))
for (i in c(1:length(intermembrane_particles$pattern))) {
  intermembrane_particles$pattern[[i]] =rescale(intermembrane_particles$pattern[[i]],512)
  tempv[i] = intermembrane_particles$pattern[[i]]$n
}

control = intermembrane_particles$pattern[intermembrane_particles$group == "control"]
acid= intermembrane_particles$pattern[intermembrane_particles$group == "acid"]
rotenone = intermembrane_particles$pattern[intermembrane_particles$group == "rotenone"]

# this function will run p.thinnig randomly until we have set of point patterns of a given size
generate_pp = function(ppp, size,prob){
  nsim_temp = table(sample(c(1:length(ppp)), size=size, replace =T))
  temp_data =rthin(ppp[[1]],P=prob,nsim = nsim_temp[1] )
  for (qq in c(2:length(nsim_temp))) {
    temp_data =  c(temp_data , rthin(ppp[[qq]],P=prob,nsim = nsim_temp[qq] ))
  }
  return(sample(temp_data))
}

control_data = generate_pp(control,size = 5000,prob = 0.7) 
acid_data = generate_pp(acid,size = 5000,prob = 0.7) 
rotenone_data = generate_pp(rotenone,size = 5000,prob = 0.7) 
saveRDS(control_data,"membrane_control_data.Rdata")
saveRDS(acid_data,"membrane_acid_data.Rdata")
saveRDS(rotenone_data,"membrane_rotenone_data.Rdata")
