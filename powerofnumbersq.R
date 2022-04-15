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

registerDoParallel(3)
squares = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))
mm=20
timeofpowertest  =c(1:6)
squares = list(c(3,1),c(2,2),c(3,2),c(4,2),c(3,3),c(3,4))
#setwd("/home/au591455/Rstuff/Results") 
setwd("C:/Users/simon/Desktop/TestR") 
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

poweroftest = function(Outlier,Data,name,m,squares,newlog = F,minpoints=20){
  n = length(Outlier)
#logpath = "/home/au591455/Rstuff/Results/log.txt"
logpath = "C:/Users/simon/Desktop/TestR/log.txt"
path_of_log<-file(logpath)
if (newlog){
  writeLines("Start UP", path_of_log)
}

  
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
  # 
  # allcomb <- is.null(nperm)
  # ncomb <- if (m1 == m2) choose(m-1, m1-1) else choose(m, m1)
  # if nperm is larger than the actual number of combinations, also use all of them
  # if (!allcomb)
  # {
  # ncomb <- if (m1 == m2) choose(m-1, m1-1) else choose(m, m1)
  #  if (ncomb < (nperm + 1)) allcomb <- TRUE
  # }
  # if (allcomb)
  #   {
  #   if (m1 == m2) index1 <- rbind(1, combn(m - 1, m1 - 1) + 1)
  #   else index1 <- combn(m, m1)
  # } else {
  #   if (m1 == m2) index1 <- rbind(1, replicate(nperm, sample(m - 1, m1 - 1) + 1)) 
  #   else index1 <- replicate(nperm, sample(m, m1)) 
  #   index1 <- cbind((1 : m1), index1) # the first is the original
  # }
  # 
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

  
  
  OutlierPPP_Permu = function(Outlier,PPP,nx,ny=nx,minpoints=20,use.tbar=FALSE,rinterval,nperm=999,sumfunc=Kest,...){
    
    grid1 = quadrats(Outlier,nx=nx,ny=ny)
    splitOutlier = split(Outlier,f=grid1)
    OutlierStat = NULL
    for(i in c(1:(nx*ny))){
      if(splitOutlier[[i]]$n >= minpoints){
        if (is.null(OutlierStat)){
          TEMPF =  sumfunc(splitOutlier[[i]],r=rinterval)
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
#### ==== Comment from Ute: where does j occur in this inner loop? =========        
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
  

  rinterval = seq(0,0.125,length.out = 30)
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %:%
    foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
      QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))],rinterval = rinterval , nx=squares[[j]][1],ny=squares[[j]][2],minpoints = minpoints)
      QQQ$p.value}
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    print(q)
    #setwd("/home/au591455/Rstuff/Results") 
    setwd("C:/Users/simon/Desktop/TestR") 
    tempC = readLines(path_of_log)
    writeLines(c(tempC,paste(name,"beginning the ",q," part ", Sys.time())),path_of_log)
    ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat")) %:%
      foreach (j= c(1:length(squares)), .combine="c", .packages = c("spatstat")) %dopar% {
        QQQ = OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*20-19):(i*20))],rinterval = rinterval, nx=squares[[j]][1],ny=squares[[j]][2],minpoints =minpoints)
        QQQ$p.value
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
tic()
poweroftest(Outlier = Matern4a,Data=Data,name ="PowerMaternA.Rdata",m=mm,squares = squares,newlog=T )
T1 = toc()
timeofpowertest[1] = T1$toc-T1$tic

Matern4b = readRDS(file = "Matern_b.Rdata")
tic()
poweroftest(Outlier = Matern4b,Data=Data,name ="PowerMaternB.Rdata",m=mm,squares = squares)
T1 = toc()
timeofpowertest[2] = T1$toc-T1$tic

Clust4a = readRDS(file = "Clust_a.Rdata")

poweroftest(Outlier = Clust4a ,Data=Data,name ="PowerClusterA.Rdata",m=mm,squares = squares )

Clust4b = readRDS(file = "Clust_b.Rdata")

poweroftest(Outlier = Clust4b ,Data=Data,name ="PowerClusterB.Rdata",m=mm,squares = squares)


Clust4c = readRDS(file = "Clust_c.Rdata")

poweroftest(Outlier = Clust4c ,Data=Data,name ="PowerClusterC.Rdata",m=mm,squares = squares )



poistest  = readRDS(file = "poisPPP.Rdata")

poweroftest(Outlier = poistest[1:50] ,Data=Data[1:1000],name ="PowerpoisTest.Rdata",m=mm,squares = squares,minpoints = 4 )

poistestP  = readRDS(file = "PowerpoisTest.Rdata")

poistestP2  = readRDS(file = "PowerpoisTest2.Rdata")

poistestP3  = readRDS(file = "PowerpoisTest3.Rdata")


rpoispp(100,nsim=nn)
tempwin =   owin(xrange=c(0,2), yrange=c(0,2))
poweroftest(Outlier = rpoispp(100,win=tempwin,nsim=50),Data=rpoispp(100,win=tempwin,nsim=50*20),name ="PowerpoisTest4.Rdata",m=mm,squares = squares,minpoints = 20 )
poistestP4  = readRDS(file = "PowerpoisTest4.Rdata")
rowMeans(poistestP4< 0.05)

tempwin =   owin(xrange=c(0,3), yrange=c(0,3))
poweroftest(Outlier = rpoispp(100,win=tempwin,nsim=100),Data=rpoispp(100,win=tempwin,nsim=100*20),name ="PowerpoisTest5.Rdata",m=mm,squares = squares,minpoints = 20 )
poistestP5  = readRDS(file = "PowerpoisTest5.Rdata")
rowMeans(poistestP5< 0.05)



poweroftest(Outlier = rpoispp(100,nsim=100),Data=rpoispp(100,nsim=100*20),name ="PowerpoisTest5.Rdata",m=mm,squares = squares,minpoints = 20 )
poistestP5  = readRDS(file = "PowerpoisTest5.Rdata")
rowMeans(poistestP5< 0.05)
mean(poistestP5[3,][-which(is.na(poistestP5[3,]))]<0.05)
mean(poistestP5[4,][-which(is.na(poistestP5[4,]))]<0.05)
mean(poistestP5[5,][-which(is.na(poistestP5[5,]))]<0.05)

poweroftest(Outlier = rpoispp(200,nsim=100),Data=rpoispp(200,nsim=100*20),name ="PowerpoisTest6.Rdata",m=mm,squares = squares,minpoints = 20 )
poistestP6  = readRDS(file = "PowerpoisTest6.Rdata")
rowMeans(poistestP6< 0.05)

OT = rpoispp(250,nsim=100)
DA= rpoispp(250,nsim=100*20)

poweroftest(Outlier = OT,Data=DA,name ="PowerpoisTest7.Rdata",m=mm,squares = squares,minpoints = 20 )
poistestP7  = readRDS(file = "PowerpoisTest7.Rdata")
rowMeans(poistestP7< 0.05)
stopImplicitCluster()


