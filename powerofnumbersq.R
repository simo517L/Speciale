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
setwd("/home/au591455/Rstuff/Results") 
#setwd("C:/Users/simon/Desktop/TestR") 
registerDoParallel(2)
nn =10000
mm=20
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

logpath = "/home/au591455/Rstuff/Results/log.txt"
fileConn<-file(logpath)
writeLines("Start UP", fileConn)
poweroftest = function(Outlier,Data,name,n,m,squares,path_of_log){
  studpermut.test.Ute <- function (foos1, foos2, use.tbar=FALSE, nperm = 25000)
  {
    ##### preparations ----------------
    n <- dim(foos1)[1]
    m1 <- dim(foos1)[2]
    stopifnot (m1 > 1) # need at least two per group
    stopifnot(dim(foos2)[1] == n) # dimensions
    m2 <- dim(foos2)[2]
    stopifnot (m2 > 1) # need at least two per group
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
  
  
  
  OutlierPPP_Permu = function(Outlier,PPP,nx,ny=nx,minpoints=20,use.tbar=FALSE,rinterval=NULL,nperm=999,sumfunc=Kest,...){
    
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
    
    return(studpermut.test.Ute(foos1 = OutlierStat,foos2= PPPStat,use.tbar=use.tbar,nperm=nperm))
  }
  
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %:%
    foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
      QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))], nx=squares[[j]][1],ny=squares[[j]][2],minpoints = 4)
      QQQ$p.value}
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    ResultpowerM1temp2 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))], nx=squares[[j]][1],ny=squares[[j]][2],minpoints =4)
        QQQ$p.value
      }
    if (!prod(ResultpowerM1temp2 == ResultpowerM1temp1)){
      ResultpowerM1  = cbind(ResultpowerM1,ResultpowerM1temp2)
      
      saveRDS(ResultpowerM1,file = name)
      ResultpowerM1temp1 = ResultpowerM1temp2
    } else{
      TEMPM = matrix(NaN,nrow = 6,ncol=5)
      ResultpowerM1  = cbind(ResultpowerM1,TEMPM)
      
      saveRDS(ResultpowerM1,file = name)
    }
    tempC = readLines(path_of_log)
    procent = q*m/n
    writeLines(c(tempC,paste(name,":",procent," % done")),path_of_log)
    }
}



tic()
poweroftest(Outlier = Matern4a,Data=Data,name ="PowerMaternA.Rdata",n=nn,m=mm,squares = squares,path_of_log = fileConn  )
T1 = toc()
timeofpowertest[1] = T1$toc-T1$tic

tic()
poweroftest(Outlier = Matern4b,Data=Data,name ="PowerMaternB.Rdata",n=nn,m=mm,squares = squares,path_of_log = fileConn  )
T1 = toc()
timeofpowertest[2] = T1$toc-T1$tic


tic()
poweroftest(Outlier = Clust4a ,Data=Data,name ="PowerClusterA.Rdata",n=nn,m=mm,squares = squares,path_of_log = fileConn  )
T1 = toc()
timeofpowertest[3] = T1$toc-T1$tic

tic()
poweroftest(Outlier = Clust4b ,Data=Data,name ="PowerClusterB.Rdata",n=nn,m=mm,squares = squares ,path_of_log = fileConn )
T1 = toc()
timeofpowertest[4] = T1$toc-T1$tic


tic()
poweroftest(Outlier = Clust4c ,Data=Data,name ="PowerClusterC.Rdata",n=nn,m=mm,squares = squares ,path_of_log = fileConn )
T1 = toc()
timeofpowertest[5] = T1$toc-T1$tic


tic()
poweroftest(Outlier = poistest ,Data=Data,name ="Powerpois.Rdata",n=nn,m=mm,squares = squares,path_of_log = fileConn  )
T1 = toc()
timeofpowertest[6] = T1$toc-T1$tic

save(timeofpowertest,file = "timeofpowertest.Rdata")
close(fileConn)
stopImplicitCluster()


