library(spatstat)
library(dplyr)
library(ppMeasures)
library(cluster)
# This is a wrapper function for stDist, so it works on point patterns. 
stDistPP = function(X,Y,...){
  p1 = cbind(X$x,X$y)
  p2 = cbind(Y$x,Y$y)
  return(stDist(p1, p2,alg="IMA", ...))
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


# Creates a counter example for the sym. nearest point distacne being a metric
set.seed(1)
windowA = owin(c(-0.5,2),c(-0.5,1))
XP = rpoint(n=7, win= disc(radius=0.1, centre=c(0,0)))
YP = ppp(x=c(1),y=c(0),window = windowA)
ZP = ppp(x=c(0),y=c(0.5),window = windowA)
marks(XP) <- factor("x")
marks(YP) <- factor("y")
marks(ZP) <- factor("z")

plot(superimpose(XP,YP,ZP), main = "")

nearest_point_metric(XP,ZP) + nearest_point_metric(ZP,YP) - nearest_point_metric(XP,YP)







# runs a loop to find a seed which has a counter example
for (i in c(1:100)){
  set.seed(i)
  X = rpoint(100)
  Y = rpoint(100)
  Z = rpoint(100)
  MM = distMppp(list(X,Y,Z),nx=2,ny=2,method=1,minpoints=20,sumfunc=Kest)
  Result =  MM[1,2] + MM[2,3] - MM[1,3]
  if(Result < 0){
    print(c("SUCCES",i))
    XRL = X
    YRL = Y
    ZRL = Z
  }
}

# Counterexample for the  triangle inequality for T(X,Y)
set.seed(3)
XRL = rpoint(100)
YRL = rpoint(100)
ZRL = rpoint(100)
MM = distMppp(list(XRL,YRL,ZRL),nx=2,ny=2,method=1,minpoints=20,sumfunc=Kest) # using our distance function, we calculate the distance matrix for the T stat.
MM[1,2] + MM[2,3] - MM[1,3] # = -0.671362
Counter_Example2 = hyperframe(list(XRL,YRL,ZRL),row.names=c("X","Y","Z")) # Collects the three point patterns in a hyper frame to make is easier to plot them.
plot(Counter_Example2,main = " ")

