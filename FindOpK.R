#We wish it fine the optimal K values for LOF when used on point patterns
library(spatstat)
library(utils)
library(ppMeasures)

# We define the function, there will be needed. 

# This is a wrapper function for stDist, so it works on point patterns. 
stDistPP = function(X,Y,...){
  p1 = cbind(X$x,X$y)
  p2 = cbind(Y$x,Y$y)
  return(stDist(p1, p2,pm=1,by,alg="IMA",bypassCheck=T, ...))
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


#This function will calculate the LOF value over a interval of $k$'s, for some distance matrix made from a set of point patterns $X$.
outlier_factors_PP = function(X,k,nx,ny=ny,method=1,minpoints=20,...){
  M = distMppp(X,nx=nx,ny=ny,method=method,minpoints=minpoints,...)
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

# We make a function, to make a similar plot like in the document
plot_k = function(X,kval,title = ""){
  max_val = do.call(pmax,as.data.frame(t(X)))
  min_val = do.call(pmin,as.data.frame(t(X)))
  mean_val = colMeans(X)
  sdev = apply(X,2,sd)
  plot(kval,mean_val,main=title,ylim = c(min(min_val),max(max_val)),xlab="k",ylab="Outlier factor (LOF)")
  arrows(kval, mean_val-sdev, kval, mean_val+sdev, length=0.01, angle=90, code=3) #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
  points(kval,max_val,  col=2 )
  points(kval,min_val,col=3)
  
  }

# We make the plot with the spike distance function
setwd("C:/Users/simon/Desktop/TestR")
Data =  readRDS(file = "DataPPP.Rdata")
OF_st = outlier_factors_PP(Data[1:100],method = 3,nx=2,ny=2,minpoints = 5,k = c(2:50))
#OF_st = outlier_factors_PP(Data[101:200],method = 3,nx=2,ny=2,minpoints = 5,k = c(2:50))

plot_k(OF_st,c(2:50))
legend("right",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_st,2,sd)
plot(sdev)



# We make the plot with the nearest neighbor function 
OF_nn = outlier_factors_PP(Data[1:100],method = 2,nx=2,ny=2,minpoints = 5,k = c(2:50))

plot_k(OF_nn, c(2:50))
legend("topright",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_nn,2,sd)
plot(sdev)

OF_T4 = outlier_factors_PP(Data[1:100],method = 1,nx=2,ny=2,minpoints = 5,k = c(2:50))
OF_T9 = outlier_factors_PP(Data[1:100],method = 1,nx=3,ny=3,minpoints = 5,k = c(2:50))
OF_T6 = outlier_factors_PP(Data[1:100],method = 1,nx=2,ny=3,minpoints = 5,k = c(2:50))
OF_T12 = outlier_factors_PP(Data[1:100],method = 1,nx=3,ny=4,minpoints = 5,k = c(2:50))


plot_k(OF_T4, c(2:50))
plot_k(OF_T6, c(2:50))
plot_k(OF_T9, c(2:50))
plot_k(OF_T12, c(2:50))
legend("topright",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_T4,2,sd)
plot(sdev)



plot_k(OF_T9, c(2:50))
legend("topright",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_T9,2,sd)
plot(sdev)

#https://stackoverflow.com/questions/10389967/common-legend-for-multiple-plots-in-r
par(oma = c(4,0,0,0), mfrow = c(2, 3), mar = c(4, 4, 1.5, 0.5))
plot_k(OF_st,c(2:50),title = "Spike-time")
#legend("right",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))

plot_k(OF_nn, c(2:50),title = "Nearest-point")
plot_k(OF_T4, c(2:50),title = "T-stat. with 4 squares")
plot_k(OF_T6, c(2:50),title = "T-stat. with 6 squares")
plot_k(OF_T9, c(2:50),title = "T-stat. with 9 squares")
plot_k(OF_T12, c(2:50),title = "T-stat. with 12 squares")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",cex =0.75, legend =c("max","mean with std","min"),pch = c(1,1,1),col=c("red","black","green"),horiz = TRUE)

par(mfrow=c(1,1))



#makes plot to show problem with small k value

TESTPPP = clickppp()
TESTDIST = crossdist(TESTPPP,TESTPPP)
k=c(2,3,4,10)
n = TESTPPP$n
m = length(k)
Result = matrix(0,nrow = n,ncol = m)
colnames(Result) <- k
for (i in  c(1:m)){
  for (j in c(1:n)){
    Result[j,i] = l_outlier_factor(TESTDIST ,k[i],j)
  }}


par(mfrow=c(1,2))
marks4 = Result[,2]
#marks4[marks4<1] = 1
marks(TESTPPP) = marks4
plot(TESTPPP,main = "k=3")
marks4 = Result[,4]
#marks4[marks4<1] = 1
marks(TESTPPP) = marks4
plot(TESTPPP,main = "k=10")
par(mfrow=c(1,1))


#makes plot to show problem with big k value

TESTPPP2 = clickppp()
TESTDIST2 = crossdist(TESTPPP2,TESTPPP2)
k=c(5,10,20)
n = TESTPPP2$n
m = length(k)
Result = matrix(0,nrow = n,ncol = m)
colnames(Result) <- k
for (i in  c(1:m)){
  for (j in c(1:n)){
    Result[j,i] = l_outlier_factor(TESTDIST2 ,k[i],j)
  }}


par(mfrow=c(1,2))
marks4 = Result[,1]

marks(TESTPPP2) = marks4
plot(TESTPPP2,main = "k=5")
marks4 = Result[,2]

marks(TESTPPP2) = marks4
plot(TESTPPP2,main = "k=10")
par(mfrow=c(1,1))

