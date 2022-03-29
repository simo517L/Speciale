
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



studpermut.test.Ute <- function (foos1, foos2, use.tbar=FALSE, nperm = 25000){
  ##### preparations ----------------
  if (is.null(foos1) |  is.null(foos2) ){
    print("foos1 or foos2 are null")
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
    print("foos1 is not a matrix")
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
    print("foos2 does not have the same lenght as foos1")
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
    print("foos2 is not a matrix")
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
  # get the permutations. 
  # If m1 == m2, break the symmetry and save half time and memory!
  
  allcomb <- is.null(nperm)
  ncomb <- if (m1 == m2) choose(m-1, m1-1) else choose(m, m1)
  # if nperm is larger than the actual number of combinations, also use all of them
  if (!allcomb)
  {
    # ncomb <- if (m1 == m2) choose(m-1, m1-1) else choose(m, m1)
    if (ncomb < (nperm + 1)) allcomb <- TRUE
  }
  if (allcomb) 
  {
    if (m1 == m2) index1 <- rbind(1, combn(m - 1, m1 - 1) + 1)
    else index1 <- combn(m, m1)
  } else {
    if (m1 == m2) index1 <- rbind(1, replicate(nperm, sample(m - 1, m1 - 1) + 1)) 
    else index1 <- replicate(nperm, sample(m, m1)) 
    index1 <- cbind((1 : m1), index1) # the first is the original
  }
  
  # do the calculations the good old fashioned way with sums and sqs, to save time
  
  SX. <- apply (foos, 1, sum)
  SXX. <- apply (foos^2, 1, sum)
  
  Tstatistic <- function (ind) # could be further optimized in symmetric case 
  {
    SX1 <- apply(foos[, ind], 1, sum)
    SXX1 <- apply(foos[, ind]^2, 1, sum)
    SX2 <- SX. - SX1
    SXX2 <- SXX. - SXX1
    mu1 <- SX1 / m1 
    mu2 <- SX2 / m2
    ss1 <- (SXX1 - (SX1^2 / m1)) / ((m1-1) / m1)
    ss2 <- (SXX2 - (SX2^2 / m2)) / ((m2-1) / m2)
    
    if (use.tbar) return (sum((mu1 -mu2)^2) / sum((ss1 + ss2))) else 
      return (mean((mu1 -mu2)^2 / (ss1 + ss2), na.rm=T))
  }
  
  Tvals <- apply(index1, 2, Tstatistic)
  
  pval <- mean(Tvals >= Tvals[1])           
  stat <- Tvals[1]
  names(stat) <- if(use.tbar) "Tbar" else "T"
  datname <- paste( deparse(substitute(foos1)),"and", deparse(substitute(foos2)))
  method <- c(paste("Studentized two sample permutation test for fda, using T",
                    ifelse(use.tbar, "bar", ""), sep=""),
              ifelse(allcomb, paste("exact test, using all",ncomb,"permutations (combinations)"), 
                     paste("using",nperm,"randomly selected permuations")))
  alternative <- "samples not exchangeable"
  ptt <- list(statistic = stat, 
              p.value = pval, 
              alternative = alternative, 
              method = method, 
              data.name = datname)
  class(ptt) <- "htest"
  return(ptt)
}

OutlierPPP_Permu = function(Outlier,PPP,nx,ny=nx,minpoints=20,use.tbar=1,nperm=1,rinterval=NULL,sumfunc=Kest,...){
  grid1 = quadrats(Outlier,nx=nx,ny=ny)
  splitOutlier = split(Outlier,f=grid1)
  OutlierStat = NULL
  for(i in c(1:(nx*ny))){
    if(splitOutlier[[i]]$n >= minpoints){
      if (is.null(OutlierStat)){
        TEMPF =  sumfunc(splitOutlier[[i]],r=rinterval)
        if (is.null(rinterval)){
          rinterval = TEMPF$r
        }
        OutlierStat = matrix(TEMPF$iso , byrow = F, ncol = 1,nrow = length(TEMPF$iso))
      } else{
        OutlierStat = cbind(OutlierStat,sumfunc(splitOutlier[[i]],r=rinterval)$iso )
      }
    }
  }
  n = length(PPP)
  PPPStat = NULL
  if (is.ppplist(PPP)){
    splitPP = list()
    for (i in c(1:n)){
      grid2 = quadrats(PPP[[i]],nx=nx,ny=ny)
      splitPP = append(splitPP,split(PPP[[i]],f=grid2))
      for(j in c(1:(nx*ny))){
        if(splitPP[[i]]$n >= minpoints){
          if (is.null( PPPStat)){
            PPPStat = matrix(sumfunc(splitPP [[i]],r=rinterval)$iso , byrow = F, ncol = 1,nrow = length(sumfunc(splitPP[[i]],r=rinterval)$iso))
          } else{
            PPPStat = cbind(PPPStat,sumfunc(splitPP [[i]],r=rinterval)$iso )
          }
        }
      }
    }
  } else{
    grid1 = quadrats(PPP,nx=nx,ny=ny)
    splitPPP = split(PPP,f=grid1)
    for(i in c(1:(nx*ny))){
      if(splitPPP[[i]]$n >= minpoints){
        if (is.null(PPPStat)){
          TEMPF =  sumfunc(splitPPP[[i]],r=rinterval)
          PPPStat= matrix(TEMPF$iso , byrow = F, ncol = 1,nrow = length(TEMPF$iso))
        } else{
          PPPStat = cbind(PPPStat,sumfunc(splitPPP[[i]],r=rinterval)$iso )
        }
      }
    }
    
    
  }
  
  
  
  return(studpermut.test.Ute(foos1 = OutlierStat,foos2= PPPStat,use.tbar=use.tbar,nperm=nperm))
}

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
nearest_point_metric = function(X,Y){
  return(nearest_pointdist(X,Y)+nearest_pointdist(Y,X))
}

distMppp = function(X,nx=3,ny=nx,method=1,minpoints=20,sumfunc=Kest){
  
  n = length(X)
  M = matrix(0,n,n)
  q=1
  for (i in c(1:n)){
    print(i)
    for (j in c(q:n)){
      if (method==1 ){
        tempstore = OutlierPPP_Permu(X[[i]],X[[j]],nx=nx,ny=ny,minpoints=minpoints,use.tbar=1,sumfunc=sumfunc,nperm=1)
        M[i,j]=tempstore$statistic
      } else if (method==2){
        M[i,j]=nearest_point_metric(X[[i]],X[[j]])
      } 
      
      M[j,i]=M[i,j]
    }
    q = q+1
  }
  return(M)
}


load("megacaryocytes.rda")

load("intermembrane_particles.rda")
plot(megacaryocytes[1:10])
plot(megacaryocytes[11:20])
intermembrane_particles
plot(intermembrane_particles[1:10])
plot(intermembrane_particles[34:44])
plot(intermembrane_particles[48:58])
TEMPM = distMppp(intermembrane_particles$pattern,method = 2)
n  =dim(TEMPM)[1]
countTRI = 0
for (i in c(1:n)){
  print(i)
    for (j in c(1:n)){
    for (q in c(1:n)){
      if (TEMPM[i,j] >  TEMPM[i,q] + TEMPM[q,j]){
        countTRI = countTRI+1}
    }
  } 
}
  
library(cluster)
hctest1= agnes(TEMPM)
plot(hctest1)

plot(agnes(TEMPM,method = "average"))
plot(agnes(TEMPM,method = "single"))
plot(agnes(TEMPM,method = "complete"))

cutree(agnes(TEMPM,method = "complete"),k=3)
TEMPM2 = distMppp(intermembrane_particles$pattern,method = 1,nx=3,ny=3)
cor1 = c(1:68)[-which(colSums(is.nan(TEMPM2) == T) == 68)]
TEMPM2a =TEMPM2[cor1,cor1]
hctest2= agnes(TEMPM2a)
plot(hctest2)
n = dim(TEMPM2a)[1]
countTRI = 0
for (i in c(1:n)){
  print(i)
  for (j in c(1:n)){
    for (q in c(1:n)){
      if (TEMPM2a[i,j] >  TEMPM2a[i,q] + TEMPM2a[q,j]){
        ix=i
        jx=j
        qx=q
        countTRI = countTRI+1}
    }
  } 
}
TEMPM2a[ix,jx] 
TEMPM2a[ix,qx] + TEMPM2a[qx,jx]

TEMPM3 = distMppp(intermembrane_particles$pattern,method = 1,nx=2,ny=2)

colSums(is.nan(TEMPM3) == T)


hctest3= agnes(TEMPM3)
plot(hctest3)


TEMPM4 = distMppp(intermembrane_particles$pattern,method = 1,nx=3,ny=3, minpoints = 10)

colSums(is.nan(TEMPM4) == T)


hctest4= agnes(TEMPM4)
plot(hctest4)
