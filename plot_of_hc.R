
library(spatstat)
library(parallel)
library(codetools )
library(iterators)
library(tictoc)
library(foreach)
library(doParallel)
library(cluster)
library(utils)

library(ppMeasures)

library(foreign)
library(maps)
library(shapefiles)
library(sp)
library(fossil)

library(colorspace)
library(dendextend)
library(circlize)
# We define the function, there will be needed. 

# This is a wrapper function for stDist, so it works on point patterns. 
stDistPP = function(X,Y,...){
  p1 = cbind(X$x,X$y)
  p2 = cbind(Y$x,Y$y)
  return(stDist(p1, p2,pm=1,alg="IMA", bypassCheck=T, ...))
}

# This is function supplied by Ute Hahn to run the permutation test for the T stat. 
studpermut.test.Ute <- function (foos1, foos2, use.tbar=FALSE, nperm = 25000){
  ##### preparations ----------------
  if (is.null(foos1) |  is.null(foos2) ){
    ptt <- list(statistic = NaN, 
                p.value = NaN, 
                alternative = "foos1 or foos2 are null", 
                method = "No method", 
                data.name = "Null")
    class(ptt) <- "htest"
    return(ptt)
  }
  n <- dim(foos1)[1]
  m1 <- dim(foos1)[2]
  if (m1 < 2){
    datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
    ptt <- list(statistic = NaN, 
                p.value = NaN, 
                alternative = "foos1 is not a matrix", 
                method = "No method", 
                data.name = datname)
    class(ptt) <- "htest"
    return(ptt)
    
    
  } # need at least two per group
  if(dim(foos2)[1] != n){
    datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
    ptt <- list(statistic = NaN, 
                p.value = NaN, 
                alternative = "foos2 does not have the same lenght as foos1", 
                method = "No method", 
                data.name = datname)
    class(ptt) <- "htest"
    return(ptt)
  } # dimensions
  m2 <- dim(foos2)[2]
  if (m2 < 2){
    datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
    ptt <- list(statistic = NaN, 
                p.value = NaN, 
                alternative = "foos2 is not a matrix", 
                method = "No method", 
                data.name = datname)
    class(ptt) <- "htest"
    return(ptt)
  } # need at least two per group
  
  m <- m1+m2
  foos <- cbind(foos1, foos2)
  index1 <- cbind(1:m1, replicate(nperm, sample(m, m1)))
  
  # do the calculations the good old fashioned way with sums and sqs, to save time
  
  foos_sq <- foos^2
  SX. <- apply (foos, 1, sum)
  SXX. <- apply (foos_sq, 1, sum)
  
  Tstatistic <- function (ind) # 
  {
    SX1 <- apply(foos[, ind], 1, sum)
    SXX1 <- apply(foos_sq[, ind], 1, sum)
    SX2 <- SX. - SX1
    SXX2 <- SXX. - SXX1
    mu1 <- SX1 / m1 
    mu2 <- SX2 / m2
    ss1 <- (SXX1 - (SX1^2 / m1)) / ((m1-1) * m1)
    ss2 <- (SXX2 - (SX2^2 / m2)) / ((m2-1) * m2)
    
    ss <- ss1 + ss2
    meandiff_sq <- (mu1 - mu2)^2
    
    result <- if (use.tbar) (sum(meandiff_sq) / sum(ss)) 
    else (mean(meandiff_sq / ss, na.rm=T))
    result
  }
  
  Tvals <- apply(index1, 2, Tstatistic)
  
  pval <- mean(Tvals >= Tvals[1], na.rm = TRUE)           
  stat <- Tvals[1]
  names(stat) <- if(use.tbar) "Tbar" else "T"
  datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
  method <- c(paste("Studentized two sample permutation test for fda, using T",
                    ifelse(use.tbar, "bar", ""), sep=""),
              paste("using",nperm,"randomly selected permuations"))
  alternative <- "samples not exchangeable"
  ptt <- list(statistic = stat, 
              p.value = pval, 
              alternative = alternative, 
              method = method, 
              data.name = datname)
  class(ptt) <- "htest"
  return(ptt)
}


# This is the function we use to calculate the T stat. distance for the distMPPP
Permu_dist = function(PPP1,PPP2,nx,ny=nx,minpoints=20,use.tbar=FALSE,rinterval,nperm=1,sumfunc=Kest,...){
  
  grid1 = quadrats(PPP1,nx=nx,ny=ny)
  splitPPP1 = split(PPP1,f=grid1)
  ResultPPP1 = NULL
  for(i in c(1:(nx*ny))){
    if(splitPPP1[[i]]$n >= minpoints){
      if (is.null(ResultPPP1 )){
        TEMPF =  sumfunc(splitPPP1[[i]],r=rinterval)
        ResultPPP1 = matrix(TEMPF$iso , byrow = F, ncol = 1,nrow = length(TEMPF$iso))
      } else{
        ResultPPP1 = cbind(ResultPPP1 ,sumfunc(splitPPP1[[i]],r=rinterval)$iso )
      }
    }
  }
  grid2 = quadrats(PPP2,nx=nx,ny=ny)
  splitPPP2 = split(PPP2,f=grid2)
  ResultPPP2 = NULL
  for(i in c(1:(nx*ny))){
    if(splitPPP2[[i]]$n >= minpoints){
      if (is.null(ResultPPP2)){
        TEMPF =  sumfunc(splitPPP2[[i]],r=rinterval)
        ResultPPP2 = matrix(TEMPF$iso , byrow = F, ncol = 1,nrow = length(TEMPF$iso))
      } else{
        ResultPPP2 = cbind(ResultPPP2,sumfunc(splitPPP2[[i]],r=rinterval)$iso )
      }
    }
  }
  
  return(studpermut.test.Ute(foos1 = ResultPPP1,foos2= ResultPPP2,use.tbar=use.tbar,nperm=nperm))
}

# This function calculate the nearest point distance.
nearest_pointdist = function(X,Y){
  LX = coords(X)
  LY = coords(Y)
  FF = function(x){
    sqrt(x[1]^2 +x[2]^2)
  }
  result = 0
  n = X$n
  for (i in c(1:n)){
    result = result+  min(apply(LY - c(LX[i,1:2]),1 , FUN = FF))
  }
  return(result)
}
#This is function used to calculate the symmetrical nearest point distance
nearest_point_metric = function(X,Y){
  return(nearest_pointdist(X,Y)+nearest_pointdist(Y,X))
}
#The distMppp function was made to calculate distance matrices for a set of point patterns
distMppp = function(X,nx=3,ny=nx,method=1,minpoints=20,sumfunc=Kest,rinterval = seq(0,0.125,length.out = 30),...){
  n = length(X)
  M = matrix(0,n,n)
  q=1
  for (i in c(1:n)){
    for (j in c(q:n)){
      if (method==1 ){
        tempstore = Permu_dist(X[[i]],X[[j]],nx=nx,ny=ny,minpoints=minpoints,use.tbar=1,sumfunc=sumfunc,rinterval = rinterval,...)
        M[i,j]=tempstore$statistic
      } else if (method==2){
        M[i,j]=nearest_point_metric(X[[i]],X[[j]],...)
      } else if (method==3){
        M[i,j]=stDistPP(X[[i]],X[[j]],...)$distance
      }
      
      M[j,i]=M[i,j]
    }
    q = q+1
  }
  return(M)
}

# We use this function to calculate the confusion matrix for two set of labels. 
conf_matrix = function(vec1,vec2){
  n = length(vec1)
  TP=0
  TN=0
  FP=0
  FN=0
  for (i in c(1:n)){
    for (j in c(1:n)){
      if (vec1[i] == vec1[j]){
        if (vec2[i] == vec2[j]){
          TP = TP+1
        } else{
          FN = FN+1
        }
      }
      else if(vec2[i] == vec2[j]){
        FP = FP+1
      } else {TN= TN +1}
      
      
      
    }
  }
  return(matrix(c(TP,FN,FP,TN),2,2,byrow = T))
}
#This function returns the evaluation value for the cluster made from a agglomerative alg.  
eval_clust = function(distmatrix,k,true_label,method="average"){
  n=length(method)
  if (n==1){
    hc = agnes(distmatrix,method = method)
    est_label = cutree(hc,k=k)
    randindex = adj.rand.index(est_label ,TrueLabel)
    CM = conf_matrix(true_label,est_label)
    Pr = CM[1,1]/(CM[1,1]+CM[2,1])  
    R = CM[1,1]/(CM[1,1]+CM[1,2])
    F_measure = 2*Pr*R/(Pr+R)
    return(c(hc$ac,randindex,F_measure)) 
  } else{
    acval = c(1:n)
    randval = c(1:n)
    Fval = c(1:n)
    for(i in c(1:n)){
      hc = agnes(distmatrix,method = method[i])
      acval[i] = hc$ac
      est_label = cutree(hc,k=k)
      randval[i]  = adj.rand.index(est_label ,TrueLabel)
      CM = conf_matrix(true_label,est_label)
      Pr = CM[1,1]/(CM[1,1]+CM[2,1])  
      R = CM[1,1]/(CM[1,1]+CM[1,2])
      Fval[i] = 2*Pr*R/(Pr+R)
    }
  }
  
  return(c(acval,randval,Fval))
}
path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path) 


# we load the data we will test on
Data =  readRDS(file = "DataPPP.Rdata")
Matern4a = readRDS(file = "Matern_a.Rdata")
Matern4b = readRDS(file = "Matern_b.Rdata")
Clust4a = readRDS(file = "Clust_a.Rdata")
Clust4b = readRDS(file = "Clust_b.Rdata")
Clust4c = readRDS(file = "Clust_c.Rdata")
poistest  = readRDS(file = "poisPPP.Rdata")

# We start by making a plot showing of the simulated point patterns 
par(mfrow = c(3,6),mar = c(0.1,0.1, 1, 0.1))
plot(Data[[1]],main = "Pois")
plot(Matern4a[[1]],main = "Matern: r=0.05")
plot(Matern4b[[1]],main = "Matern: r=0.02")
plot(Clust4a[[1]],main = "Clust: scale=0.1,mu=1")
plot(Clust4b[[1]],main = "Clust: scale=0.1,mu=4")
plot(Clust4c[[1]],main = "Clust: scale=0.5,mu=4")
plot(Data[[2]],main = "")
plot(Matern4a[[2]],main = "")
plot(Matern4b[[2]],main = "")
plot(Clust4a[[2]],main = "")
plot(Clust4b[[2]],main = "")
plot(Clust4c[[2]],main = "")
plot(Data[[3]],main = "")
plot(Matern4a[[3]],main = "")
plot(Matern4b[[3]],main = "")
plot(Clust4a[[3]],main = "")
plot(Clust4b[[3]],main = "")
plot(Clust4c[[3]],main = "")
par(mfrow = c(1,1))

# we make a function, whihc makes the dendrograms look better 
make_dend = function(hc,k,label){
  n = length( hc$order)
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, 1:n)
  dend <- color_branches(dend, k=k)
  
 # labels_colors(dend) <-
#    rainbow_hcl(k)[sort_levels_values(
#      as.numeric(as.factor(label))[order.dendrogram(dend)]
#    )]
  
  labels(dend) <- paste(as.character(label)[order.dendrogram(dend)],
                        "(",labels(dend),")", 
                        sep = "")
  #dend <- hang.dendrogram(dend,hang_height=0.1)
  dend <- set(dend, "labels_cex", 0.5)
  return(dend)
}

# Makes a example of such a dendrogram
vec = c(1:5)
MM  = distMppp(c(Data[vec],Matern4a[vec],Clust4c[vec]),method=1,nx=3,ny=3,minpoints=5)
hctest= agnes(MM)
dend = make_dend(hctest,k=3,label=c(rep("Pois",5),rep("Matern",5),rep("Cluste",5)))

plot(dend,horiz =  TRUE, nodePar = list(cex = .000001))



# Makes a dendrogram to showcase the "hard" data set
MM  = distMppp(c(Data[vec],Matern4b[vec],Clust4a[vec]),method=1,nx=2,ny=3,minpoints=5)
hctest1= agnes(MM,method = "complete")

dend1 = make_dend(hctest1,k=3,label=c(rep("Pois",5),rep("Matern",5),rep("Cluste",5)))

plot(dend1,horiz =  TRUE, nodePar = list(cex = .0001))



# we make a dendrogram to showcase the "big" dataset
vec = c(1:3)
datavec = c(Data[vec],Matern4a[vec],Matern4b[vec],Clust4a[vec],Clust4b[vec],Clust4c[vec])
MM  = distMppp(datavec,method=1,nx=2,ny=3,minpoints=5)
hctest4= agnes(MM,method = "complete")

Label = c(rep("Pois",3),rep("MaternA",3),rep("MaternB",3),rep("ClusteA",3),rep("ClusteB",3),rep("ClusteC",3))
dend4 = make_dend(hctest4,k=6,label=Label )

plot(dend4,horiz =  TRUE, nodePar = list(cex = .0001))




load(file = "intermembrane_particles.rda")
# we make the big circle dendrogram
MMfuld  = distMppp(intermembrane_particles$pattern,method=1,nx=3,ny=3,minpoints=5)
hcmembranefull= agnes(MMfuld,method = "complete")
fullLabel  = intermembrane_particles$group
dend_membranefull = as.dendrogram(hcmembranefull)

n = length( hc$order)

dend_membranefull <- color_branches(dend_membranefull, k=3)


labels(dend_membranefull) <- paste(as.character(fullLabel)[order.dendrogram(dend_membranefull)],
                      "(",labels(dend_membranefull),")", 
                      sep = "")

dend_membranefull <- hang.dendrogram(dend_membranefull,hang_height=3)
dend_membranefull <- set(dend_membranefull, "labels_cex", 0.6)
circlize_dendrogram(dend_membranefull,
                    labels_track_height = 0.2,
                    dend_track_height = 0.5)
