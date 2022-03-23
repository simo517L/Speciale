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
  if (!is.ppp(PPP)){
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
    for (j in c(q:n)){
      if (method==1 ){
        tempstore = OutlierPPP_Permu(X[[i]],X[[j]],nx=nx,ny=ny,minpoints=minpoints,use.tbar=1,sumfunc=sumfunc)
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
  
  colum = which.max(Result[n+1,])
  return(mean(Result[n+1,colum] <= Result[,colum]))
}

outlier_factors_PP(c(Matern4a[1],Data[1:20]),method = 1,nx=2,ny=2,minpoints = 4,k = c(5:15))
OutlierPPP_Permu(Matern4a[[1]],PPP=Data[2], nx = 3, ny = 3, minpoints = 10)

registerDoParallel(2)
nn =10
mm=4
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


#setwd("/home/au591455/Rstuff/Results") 
setwd("C:/Users/simon/Desktop/TestR")
powertest_OF_sq = function(Outlier,Data,name,n,m,squares,newlog = F){
  
  
  
  
  
  
  
  
  #setwd("/home/au591455/Rstuff/Results") 
  setwd("C:/Users/simon/Desktop/TestR")
  #logpath = "/home/au591455/Rstuff/Results/logOF.txt"
  logpath = "C:/Users/simon/Desktop/TestR/logOF.txt"
  fileConn<-file(logpath)
  writeLines("Start UP", fileConn)
  path_of_log = fileConn
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %:%
    foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
      Test_outlier_OF(Outlier = Outlier[[i]],PPP = Data[c((i*20-19):(i*20))],nx=squares[[j]][1],ny=squares[[j]][2],method=1,minpoints=4,Kinterval = c(10:15))
      }
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        Test_outlier_OF(Outlier = Outlier[[i]],PPP = Data[c((i*20-19):(i*20))],nx=squares[[j]][1],ny=squares[[j]][2],method=1,minpoints=4,Kinterval = c(10:15))
      }
    ResultpowerM1  = cbind(ResultpowerM1,ResultpowerM1temp1)
    saveRDS(ResultpowerM1,file = name)
    tempC = readLines(path_of_log)
    procent = q*m/n
    writeLines(c(tempC,paste(name,":",procent," % done", Sys.time())),path_of_log)
  }
  
}


Matern4a = readRDS(file = "Matern_a.Rdata")
tic()
powertest_OF_sq(Outlier = Matern4a,Data=Data,name ="PowerOFMaternA.Rdata",m=mm,squares = squares)
T1 = toc()
timeofpowertest[2] = T1$toc-T1$tic

Matern4b = readRDS(file = "Matern_b.Rdata")
tic()
powertest_OF_sq(Outlier = Matern4b,Data=Data,name ="PowerOFMaternB.Rdata",m=mm,squares = squares)
T1 = toc()
timeofpowertest[2] = T1$toc-T1$tic

Clust4a = readRDS(file = "Clust_a.Rdata")
tic()
powertest_OF_sq(Outlier = Clust4a ,Data=Data,name ="PowerOFClusterA.Rdata",m=mm,squares = squares )
T1 = toc()
timeofpowertest[3] = T1$toc-T1$tic

Clust4b = readRDS(file = "Clust_b.Rdata")
tic()
powertest_OF_sq(Outlier = Clust4b ,Data=Data,name ="PowerOFClusterB.Rdata",m=mm,squares = squares)
T1 = toc()
timeofpowertest[4] = T1$toc-T1$tic

Clust4c = readRDS(file = "Clust_c.Rdata")
tic()
powertest_OF_sq(Outlier = Clust4c ,Data=Data,name ="PowerOFClusterC.Rdata",m=mm,squares = squares )
T1 = toc()
timeofpowertest[5] = T1$toc-T1$tic

poistest  = readRDS(file = "poisPPP.Rdata")
tic()
powertest_OF_sq(Outlier = poistest ,Data=Data,name ="PowerOFpois.Rdata",m=mm,squares = squares )
T1 = toc()
timeofpowertest[6] = T1$toc-T1$tic

save(timeofpowertest,file = "timeofpowertest.Rdata")


stopImplicitCluster()
