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
library(utils)
library("ppMeasures",lib.loc=liblocation )

registerDoParallel(4)
mm=10
squares = list(c(2,2),c(3,2),c(3,3),c(3,4))

DataSize= c(20,25,30,35,40)

path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path) 
powertest_OF = function(Outlier,Data,name,n=NULL,m,squares,newlog = F,DataSize,method=1,Kinterval = c(5:10),path){
  if(min(DataSize) < Kinterval[length(Kinterval)] ){
    Kinterval = Kinterval[Kinterval< min(DataSize)]
    print("Kinterval values cannot be bigger then the amount of point patterns")
  }
  if (is.null(n)){
    n = length(Outlier)
  }
  
  #logpath = "/home/au591455/Rstuff/Results/logOF_S.txt"
  logpath = "C:/Users/simon/Desktop/SpecialeProject/Speciale/logOF_S.txt"
  path_of_log<-file(logpath)
  if (newlog){
    writeLines(paste("Start UP",Sys.time()), path_of_log)
  }

  
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
  
  # This calculates the K-distance
  K_dist = function(M,k,p){
    dist_to_point = M[p,]
    dist_to_point = sort(dist_to_point)
    return(dist_to_point[k])
  }
  # This calculates the neighborhood of a K-distance
  K_dist_neighborhood = function(M,k,p){
    K_dist_Result=K_dist(M,k,p)
    dist_to_point = M[p,]
    return(which(dist_to_point <= K_dist_Result))
  }
  # This calculates local reachability density.
  lrd = function(M,k,p){
    PP_in_neighborhood = K_dist_neighborhood(M,k,p)
    Result= 0
    for (i in c(1:length(PP_in_neighborhood))){
      Result = Result + max(K_dist(M=M,k=k,p=PP_in_neighborhood[i]),M[p,PP_in_neighborhood[i]])
    }
    return(1/( Result/length(PP_in_neighborhood))) 
  }
  # This function calculate the local outlier factor. 
  l_outlier_factor = function(M,k,p){
    PP_in_neighborhood = K_dist_neighborhood(M,k,p)
    Result= 0
    p_lrd=lrd(M,k,p) 
    for (i in c(1:length(PP_in_neighborhood))){
      Result = Result + lrd(M,k,PP_in_neighborhood[i]) / p_lrd
    }
    return( Result/length(PP_in_neighborhood)) 
  }
  # This function will calculate the LOF value over a interval of $k$'s
  outlier_factors_PP = function(X,k,nx,ny=ny,method=1,minpoints=20){
    M = distMppp(X,nx=nx,ny=ny,method=method,minpoints=minpoints)
    n = length(X)
    m = length(k)
    if (sum(is.nan(M))>0){
      Mtemp = M[-n,-n]
      while (sum(is.nan(Mtemp))>0) {
        nan_index = which.max(colSums(is.nan(Mtemp)))
        X = X[-nan_index]
        Mtemp = Mtemp[-nan_index,-nan_index]
        M = M[-nan_index,-nan_index]
      }
      if (sum(is.nan(M))>0){
        return(matrix(NaN,nrow = n,ncol = m))
      }
    }
    n = length(X)
    Result = matrix(0,nrow = n,ncol = m)
    colnames(Result) <- k
    for (i in  c(1:m)){
      for (j in c(1:n)){
        Result[j,i] = l_outlier_factor(M,k[i],j)
      }}
    return(Result)
  }
  #This function uses the LOF values to calculate a p-value. 
  Test_outlier_OF = function(Outlier,PPP,nx,ny=nx,method=1,minpoints=20,Kinterval){
    X = c(PPP,list(Outlier))
    if (length(X) < max(Kinterval)){
      Kinterval = Kinterval[Kinterval< length(X)]
      print("Kinterval values cannot be bigger then the amount of point patterns")
    }
    n = length(PPP)
    Result = outlier_factors_PP(X=X,k=Kinterval,nx=nx,ny=ny,method=method,minpoints=minpoints)
    
    Result = apply(Result,1,max)
    return(mean(Result[n+1] <= Result))
  }
  
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat","ppMeasures")) %:%
    foreach (j= c(1:length(DataSize)), .combine="c", .packages = c("spatstat","ppMeasures")) %dopar% {
      DS1 = DataSize[j]
      DS2 = DS1 -1 
      Test_outlier_OF(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))],method = method, nx=squares[1],ny=squares[2],minpoints = 5,Kinterval=Kinterval)
    }
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    setwd(path) 
    tempC = readLines(path_of_log)
    writeLines(c(tempC,paste(name,"beginning the ",q," part ", Sys.time())),path_of_log)
    ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat","ppMeasures"))%:%
      foreach (j= c(1:length(DataSize)), .combine="c", .packages = c("spatstat","ppMeasures"))  %dopar% {
        DS1 = DataSize[j]
        DS2 = DS1 -1 
        Test_outlier_OF(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))],method = method, nx=squares[1],ny=squares[2],minpoints = 5,Kinterval=Kinterval)
    }
    ResultpowerM1  = cbind(ResultpowerM1,ResultpowerM1temp1)
    saveRDS(ResultpowerM1,file = name)
    tempC = readLines(path_of_log)
    procent = q*m/n
    writeLines(c(tempC,paste(name,":",procent," % done", Sys.time())),path_of_log)
  }
  close(path_of_log)
}
Data =  readRDS(file = "DataPPP.Rdata")

Matern4a = readRDS(file = "Matern_a.Rdata")
powertest_OF(Outlier = Matern4a,Data=Data,name ="Size_MaternA_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize,newlog=T)

Matern4b = readRDS(file = "Matern_b.Rdata")
powertest_OF(Outlier = Matern4b,Data=Data,name ="Size_MaternB_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize,,newlog=T)

Clust4a = readRDS(file = "Clust_a.Rdata")
powertest_OF(Outlier = Clust4a ,Data=Data,name ="Size_ClusterA_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize, )

Clust4b = readRDS(file = "Clust_b.Rdata")
powertest_OF(Outlier = Clust4b ,Data=Data,name ="Size_ClusterB_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize,)

Clust4c = readRDS(file = "Clust_c.Rdata")
powertest_OF(Outlier = Clust4c ,Data=Data,name ="Size_ClusterC_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize, )


poistest  = readRDS(file = "poisPPP.Rdata")
powertest_OF(Outlier = poistest ,Data=Data,name ="Size_pois_OFT.Rdata",n=1000,method = 1,m=mm,squares =c(2,3),path=path,DataSize=DataSize, )




stopImplicitCluster()
