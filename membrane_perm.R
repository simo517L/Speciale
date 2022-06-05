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
#library("ppMeasures",lib.loc=liblocation )

registerDoParallel(4)
mm=20


path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path) 
powertest_membrane = function(Outlier,Data,name,n=NULL,m,squares,newlog = F,DataSize,path){
  if (is.null(n)){
    n = length(Outlier)
  }
  
  logpath = "C:/Users/simon/Desktop/SpecialeProject/Speciale/permLOG.txt"
  path_of_log<-file(logpath)
  if (newlog){
    writeLines(paste("Start UP",Sys.time()), path_of_log)
  }
  DS1 = DataSize
  DS2 = DS1 -1 
  
  
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
  
  
  
  #This function is used to run the permutation outlier test.
  OutlierPPP_Permu = function(Outlier,PPP,nx,ny=nx,minpoints=20,use.tbar=T,rinterval,nperm=999,sumfunc=Kest,...){
    
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
      splitPP = split(PPP[[i]],f=grid2)
      for(j in c(1:(nx*ny))){
        if(splitPP[[j]]$n >= minpoints){
          sumfuncR = sumfunc(splitPP[[j]],r=rinterval)$iso
          if (is.null( PPPStat)){
            PPPStat = matrix(sumfuncR , byrow = F, ncol = 1,nrow = length(sumfuncR))
          } else{
            PPPStat = cbind(PPPStat,sumfuncR )
          }
        }
      }
    }
    
    return(studpermut.test.Ute(foos1 = OutlierStat,foos2= PPPStat,use.tbar=use.tbar,nperm=nperm))
  }
  

  
  # We define the interval for the summary function
  rinterval = seq(0,0.125,length.out = 30)
  # we run the first m times, so we have some data to work with. After we keep saving after the code has run m times. 
  ResultpowerM1temp1 <- foreach (i= c(1:m), .combine="cbind", .packages = c("spatstat")) %dopar% {
    
    temp =OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))], nx=squares[1],ny=squares[2],rinterval=rinterval,minpoints=5,use.tbar=TRUE)
    temp$p.value
    }
  
  ResultpowerM1 = ResultpowerM1temp1
  for (q in c(2:(n/m))){
    setwd(path)
    tempC = readLines(path_of_log)
    writeLines(c(tempC,paste(name,"beginning the ",q," part ", Sys.time())),path_of_log)
    ResultpowerM1temp1 <- foreach (i= c((q*m-(m-1)):(q*m)), .combine="cbind", .packages = c("spatstat"))  %dopar% {
      temp =OutlierPPP_Permu(Outlier = Outlier[[i]], PPP = Data[c((i*DS1-DS2):(i*DS1))], nx=squares[1],ny=squares[2],rinterval=rinterval,minpoints=5,use.tbar=TRUE)
      temp$p.value
      }
    ResultpowerM1  = cbind(ResultpowerM1,ResultpowerM1temp1)
    saveRDS(ResultpowerM1,file = name)
    tempC = readLines(path_of_log)
    procent = q*m/n
    writeLines(c(tempC,paste(name,":",procent," % done", Sys.time())),path_of_log)
  }
  close(path_of_log)
}



control_data = readRDS("control_data.Rdata")

acid_data = readRDS("acid_data.Rdata")

rotenone_data = readRDS("rotenone_data.Rdata")

control_pop = readRDS("control_pop.Rdata")

powertest_membrane(Outlier = control_data,Data=control_pop,name ="memb_perm_control.Rdata",n=1000,m=mm,squares =c(2,3),DataSize=20,newlog=T,path=path)
powertest_membrane(Outlier = acid_data,Data=control_pop,name ="memb_perm_acid.Rdata",n=1000,m=mm,squares =c(2,3),DataSize=20,path=path)
powertest_membrane(Outlier = rotenone_data,Data=control_pop,name ="memb_perm_rotenone.Rdata",n=1000,m=mm,squares =c(2,3),DataSize=20,path=path)



stopImplicitCluster()