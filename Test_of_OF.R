liblocation = "/home/au591455/Rstuff/library"
#liblocation = NULL
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
#library("ppMeasures",lib.loc=liblocation )

registerDoParallel(4)
mm=10
squares = list(c(2,2),c(3,2),c(3,3),c(3,4))


setwd("/home/au591455/Rstuff/Results") 
#setwd("C:/Users/simon/Desktop/TestR")
powertest_OF_sq = function(Outlier,Data,name,n=NULL,m,squares,newlog = F,DataSize,method=1,Kinterval = c(5:10)){
  if(DataSize < Kinterval[length(Kinterval)] ){
  Kinterval = Kinterval[Kinterval< DataSize]
  print("Kinterval values cannot be bigger then the amount of point patterns")
}
  if (is.null(n)){
    n = length(Outlier)
  }
    
    logpath = "/home/au591455/Rstuff/Results/logOF.txt"
    #logpath = "C:/Users/simon/Desktop/TestR/logOF.txt"
    path_of_log<-file(logpath)
    if (newlog){
      writeLines(paste("Start UP",Sys.time()), path_of_log)
    }
    DS1 = DataSize
    DS2 = DS1 -1 
    
    # stDistPP = function(X,Y,...){
    #    p1 = cbind(X$x,X$y)
    #    p2 = cbind(Y$x,Y$y)
    #    return(stDist(p1, p2,alg="IMA", ...))
    #  }
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
    
    OutlierPPP_Permu =function(PPP1,PPP2,nx,ny=nx,minpoints=20,use.tbar=FALSE,rinterval,nperm=1,sumfunc=Kest,...){
      
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
    
    distMppp = function(X,nx=3,ny=nx,method=1,minpoints=5,sumfunc=Kest){
      
      n = length(X)
      M = matrix(0,n,n)
      q=1
      for (i in c(1:n)){
        for (j in c(q:n)){
          if (method==1 ){
            tempstore = OutlierPPP_Permu(X[[i]],X[[j]],nx=nx,ny=ny,minpoints=minpoints,use.tbar=1,sumfunc=sumfunc,nperm=1,rinterval = seq(0,0.125,length.out = 30))
            M[i,j]=tempstore$statistic
          } else if (method==2){
            M[i,j]=nearest_point_metric(X[[i]],X[[j]])
          } #else if (method==3){
          #  M[i,j]=nearest_point_metric(X[[i]],X[[j]])
          #}
          
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

    ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        Test_outlier_OF(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))],method = method, nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5,Kinterval=Kinterval)
        }
    
    ResultpowerM1 = ResultpowerM1temp1
    for (q in c(2:(n/m))){
      setwd("/home/au591455/Rstuff/Results") 
      #setwd("C:/Users/simon/Desktop/TestR") 
      tempC = readLines(path_of_log)
      writeLines(c(tempC,paste(name,"beginning the ",q," part ", Sys.time())),path_of_log)
      ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat")) %:%
        foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
          Test_outlier_OF(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))],method = method, nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 5,Kinterval=Kinterval)
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
#powertest_OF_sq(Outlier = Matern4a,Data=Data,name ="Power_MaternA_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20,newlog=T)

Matern4b = readRDS(file = "Matern_b.Rdata")
#powertest_OF_sq(Outlier = Matern4b,Data=Data,name ="Power_MaternB_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20)

Clust4a = readRDS(file = "Clust_a.Rdata")
#powertest_OF_sq(Outlier = Clust4a ,Data=Data,name ="Power_ClusterA_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20 )

Clust4b = readRDS(file = "Clust_b.Rdata")
powertest_OF_sq(Outlier = Clust4b ,Data=Data,name ="Power_ClusterB_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20)

Clust4c = readRDS(file = "Clust_c.Rdata")
powertest_OF_sq(Outlier = Clust4c ,Data=Data,name ="Power_ClusterC_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20 )


poistest  = readRDS(file = "poisPPP.Rdata")
powertest_OF_sq(Outlier = poistest ,Data=Data,name ="pois_OFsq.Rdata",n=1000,method = 1,m=mm,squares = squares,DataSize=20 )

stopImplicitCluster()
