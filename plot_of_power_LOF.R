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



library(xtable)
# This function uses prop.test to calcualte conf. intervals
conf_of_power=function(M,n){
  Result = rowSums(M)
  nn = length( Result)
  p_value = c(1:nn)
  upper_lim = c(1:nn)
  lower_lim = c(1:nn)
  for (i in c(1:nn)){
    TEMP = prop.test(Result[i],n)
    p_value[i] = TEMP$estimate
    upper_lim[i] = TEMP$conf.int[1]
    lower_lim[i] =  TEMP$conf.int[2]
  }
  return(data.frame(p_value=p_value,upper_lim=upper_lim,lower_lim=lower_lim))
}


# this function plots the data 
plot_of_power = function(Data,Size,X,conf = FALSE,title="",legendC = FALSE,legendNames=NULL,yval=NULL,xlab="Number of squares"){
  if (!is.null(yval)){
    a=yval[1]
    b=yval[2]
  }
  if ( is.list(Data) & length(Data) > 1){
    n = length(Data)
    Rtemp=c(1:n)
    Rlow=c(1:n)
    Rup=c(1:n)
    for (i in c(1:n)){
      Results = conf_of_power(Data[[i]],Size)
      Rtemp[i] = Results$p_value
      Rlow[i] = Results$lower_lim
      Rup[i] = Results$upper_lim
    }
    if (conf == FALSE){
      
      if (is.null(yval)){
        a = min(Rtemp) -0.05
        b = max(Rtemp) +0.05
      }
      
      plot(X,Rtemp[1],type="p",lty=1,col=1,ylim = c(a,b),xlab=xlab,ylab = "Rate of rejection",main = title)
      for ( j in c(2:n)){
        points(X,Rtemp[j],lty=j,col=j)
      }
    }     else{
      if (is.null(yval)){
        a = min(Rlow) -0.05
        b = max(Rup) +0.05
      }
      plot(X,Rtemp[1],type="p",lty=1,col=1,ylim = c(a,b),xlab=xlab,ylab = "Rate of rejection",main = title)
      arrows(X, Rlow[1],X, Rup[1], length=0.01, angle=90, code=3) 
      #points(X,Rup[1],lty=2,col=1)
      #points(X,Rlow[1],lty=2,col=1)
      for ( j in c(2:n)){
        points(X,Rtemp[j],lty=1,col=j)
        arrows(X, Rlow[j],X, Rup[j], length=0.01, angle=90, code=3)
        #points(X,Rup[j],lty=2,col=j)
        #points(X,Rlow[j],lty=2,col=j)
      }
      
    }} else if ( is.matrix(Data) & length(Data) > 1) {
      n = dim(Data)[2]
      Results = conf_of_power(Data,Size)
      
      if (conf == FALSE){
        if (is.null(yval)){
          a = min(Results$p_value) -0.1
          b = max(Results$p_value) +0.1
        }
        plot(X,Results$p_value,type="p",lty=1,col=1,xlab=xlab,ylim = c(a,b),ylab = "Rate of rejection",main = title)
        
      }     else{
        if (is.null(yval)){
          a = min(Results$lower_lim) -0.1
          b = max(Results$upper_lim) +0.1
        }
        plot(X,Results$p_value,type="p",lty=1,col=1,xlab=xlab,ylim = c(a,b),ylab = "Rate of rejection",main = title)
        arrows(X, Results$lower_lim,X, Results$upper_lim, length=0.01, angle=90, code=3)
        #points(X,Results$upper_lim,lty=2,col=1)
        #points(X,Results$lower_lim,lty=2,col=1)
        
      }} else {
        n=1
        Results = conf_of_power(Data,Size)
        if (conf == FALSE){
          if (is.null(yval)){
            a = min(Results$p_value) -0.1
            b = max(Results$p_value) +0.1
          }
          plot(X,Results$p_value,type="p",lty=1,col=1,xlab=xlab, ylim= c(a,b),ylab = "Rate of rejection",main = title)}
        else {
          if (is.null(yval)){
            a = min(Results$lower_lim) -0.1
            b = max(Results$upper_lim) +0.1
          }
          plot(X,Results$p_value,type="l",lty=1,col=1,xlab=xlab,ylim = c(a,b),ylab = "Rate of rejection",main = title)
          arrows(X, Results$lower_lim,X, Results$upper_lim, length=0.01, angle=90, code=3)
          #points(X,Results$upper_lim,lty=2,col=1)
          #points(X,Results$lower_lim,lty=2,col=1)
        }
      }
  if (legendC & n>1){
    if (is.null(legendNames)){#https://www.statology.org/legend-outside-plot-r/
      if (is.null(names(Data))){
        legend("topright" , legend = c(1:n),col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))
      } else{
        legend("topright" , legend = names(Data),col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))
      }
    } else{legend("topright" , legend = legendNames,col = 1:n,inset=c(-0.2, 0) ,pch=c(1,3))}
  }
}

path = "C:/Users/simon/Desktop/SpecialeProject/Speciale/Data" # Remember to set your own path 
setwd(path)  

# We load the results from the LOF test using the T stat. 

result_MaternA_OFsq= readRDS(file = "Power_MaternA_OFsq.Rdata")
result_MaternB_OFsq= readRDS(file = "Power_MaternB_OFsq.Rdata")
result_ClusterA_OFsq= readRDS(file = "Power_ClusterA_OFsq.Rdata")
result_ClusterB_OFsq= readRDS(file = "Power_ClusterB_OFsq.Rdata")
result_ClusterC_OFsq= readRDS(file = "Power_ClusterC_OFsq.Rdata")
result_pois_OFsq= readRDS(file = "pois_OFsq.Rdata")

number_of_squares = c(4,6,9,12) 

#Then we plot them.  
par(mfrow = c(3,2))

power_MaternA_OFsq= result_MaternA_OFsq <=  0.05
plot_of_power(Data=power_MaternA_OFsq,yval=c(0,1),title = "Hardcore II: r=0.05",Size = dim(power_MaternA_OFsq)[2],X= number_of_squares,conf = T,legendC = F)


power_MaternB_OFsq = result_MaternB_OFsq <=  0.05
sum(is.nan(colMeans(power_MaternB_OFsq)) == T)
plot_of_power(Data=power_MaternB_OFsq,yval=c(0,1),title = "Hardcore II: r=0.02",Size = dim(power_MaternB_OFsq)[2],X= number_of_squares,conf = T,legendC = F)


power_ClusterA_OFsq= result_ClusterA_OFsq <=  0.05
sum(is.nan(colMeans(power_ClusterA_OFsq)) == T)
plot_of_power(Data=power_ClusterA_OFsq,yval=c(0,1),title = "Clust: scale=0.1,mu=1",Size = dim(power_ClusterA_OFsq)[2],X= number_of_squares,conf = T,legendC = F)


power_ClusterB_OFsq= result_ClusterB_OFsq <=  0.05
sum(is.nan(colMeans(power_ClusterB_OFsq)) == T)
plot_of_power(Data=power_ClusterB_OFsq,yval=c(0,1),title = "Clust: scale=0.1,mu=4",Size = dim(power_ClusterB_OFsq)[2],X= number_of_squares,conf = T,legendC = F)



power_ClusterC_OFsq = result_ClusterC_OFsq <=  0.05

plot_of_power(Data=power_ClusterC_OFsq,yval=c(0,1),title = "Clust: scale=0.5,mu=4",Size = dim(power_ClusterC_OFsq)[2],X= number_of_squares,conf = T,legendC = F)



power_pois_OFsq = result_pois_OFsq <=  0.05
plot_of_power(Data=power_pois_OFsq,yval=c(0,1),title = "Poisson",Size = dim(power_pois_OFsq)[2],X= number_of_squares,conf = T,legendC = F)

par(mfrow = c(1,1))



# We load the results from the LOF test using the nearest point distance. 
result_MaternA_OFNP= readRDS(file = "Power_MaternA_OFNP.Rdata")
result_MaternB_OFNP= readRDS(file = "Power_MaternB_OFNP.Rdata")
result_ClusterA_OFNP= readRDS(file = "Power_ClusterA_OFNP.Rdata")
result_ClusterB_OFNP= readRDS(file = "Power_ClusterB_OFNP.Rdata")
result_ClusterC_OFNP= readRDS(file = "Power_ClusterC_OFNP.Rdata")
result_pois_OFNP= readRDS(file = "pois_OFNP.Rdata")

power_MaternA_OFNP= result_MaternA_OFNP <=  0.05
power_MaternB_OFNP = result_MaternB_OFNP <=  0.05
power_ClusterA_OFNP= result_ClusterA_OFNP <=  0.05
power_ClusterB_OFNP= result_ClusterB_OFNP <=  0.05
power_ClusterC_OFNP = result_ClusterC_OFNP <=  0.05
power_pois_OFNP = result_pois_OFNP <=  0.05

# We load the results from the LOF test using the spike-time distance. 
result_MaternA_OFST= readRDS(file = "Power_MaternA_OFST.Rdata")
result_MaternB_OFST= readRDS(file = "Power_MaternB_OFST.Rdata")
result_ClusterA_OFST= readRDS(file = "Power_ClusterA_OFST.Rdata")
result_ClusterB_OFST= readRDS(file = "Power_ClusterB_OFST.Rdata")
result_ClusterC_OFST= readRDS(file = "Power_ClusterC_OFST.Rdata")
result_pois_OFST= readRDS(file = "pois_OFST.Rdata")

power_MaternA_OFST= result_MaternA_OFST <=  0.05
power_MaternB_OFST= result_MaternB_OFST <=  0.05
power_ClusterA_OFST= result_ClusterA_OFST <=  0.05
power_ClusterB_OFST= result_ClusterB_OFST <=  0.05
power_ClusterC_OFST= result_ClusterC_OFST <=  0.05
power_pois_OFST= result_pois_OFNP <=  0.05


msd = function(x){
  temp = prop.test(x=sum(x),n=length(x))
  return(c(temp$estimate))}


rowmsd = function(x){
  t(apply(x,1,msd))
}

# we combine the data into one matrix and then create a latex table out form it
MM = rowmsd(power_MaternA_OFsq)

MM = rbind(MM,rowmsd(power_MaternB_OFsq))

MM = rbind(MM,rowmsd(power_ClusterA_OFsq))
MM = rbind(MM,rowmsd(power_ClusterB_OFsq))
MM = rbind(MM,rowmsd(power_ClusterC_OFsq))
MM = rbind(MM,rowmsd(power_pois_OFsq))

MM = cbind(MM,c(msd(power_MaternA_OFST),msd(power_MaternB_OFST),msd(power_ClusterA_OFST),
  msd(power_ClusterB_OFST),msd(power_ClusterC_OFST),msd(power_pois_OFST)))

MM = cbind(MM,c(msd(power_MaternA_OFNP),msd(power_MaternB_OFNP),msd(power_ClusterA_OFNP),
                msd(power_ClusterB_OFNP),msd(power_ClusterC_OFNP),msd(power_pois_OFNP)))

MM = data.frame(MM,row.names = c("Hardcore II: r=0.05", "Hardcore II: r=0.02", "Cluster:scale=0.1,mu=1", "Cluster:scale=0.1,mu=4", "Cluster:scale=0.5,mu=4", "Poisson"))
names(MM) = c("T (sq=4)","T (sq=6)","T (sq=9)","T (sq=12)","Spike-time","Nearest Point" )

tableMM = xtable(MM,auto=T,floating = FALSE)
print(tableMM)

# Size
# We load the result for the permutation test used on different pop. sizes. 
size_MaternA_perm= readRDS(file = "PowerMaternA_Size.Rdata") <=  0.05
size_MaternB_perm= readRDS(file = "PowerMaternB_Size.Rdata") <=  0.05
size_ClusterA_perm= readRDS(file = "PowerClusterA_Size.Rdata") <=  0.05
size_ClusterB_perm= readRDS(file = "PowerClusterB_Size.Rdata") <=  0.05
size_ClusterC_perm= readRDS(file = "PowerClusterC_Size.Rdata") <=  0.05
size_pois_perm= readRDS(file = "PowerpoisTest_Size.Rdata") <=  0.05


#We then plot them
par(mfrow = c(3,2))
DataSize= c(5,10,15,20,25,30,35,40)
plot_of_power(Data=size_MaternA_perm,yval=c(0,1),xlab = "Pop. Size",title = "Hardcore II: r=0.05",Size = dim(size_MaternA_perm)[2],X= DataSize,conf = T,legendC = F)

plot_of_power(Data=size_MaternB_perm,yval=c(0,1),xlab = "Pop. Size",title = "Hardcore II: r=0.02",Size = dim(size_MaternB_perm)[2],X= DataSize,conf = T,legendC = F)

plot_of_power(Data=size_ClusterA_perm,yval=c(0,1),xlab = "Pop. Size",title = "Clust: scale=0.1,mu=1",Size = dim(size_ClusterA_perm)[2],X= DataSize,conf = T,legendC = F)

plot_of_power(Data=size_ClusterB_perm,yval=c(0,1),xlab = "Pop. Size",title = "Clust: scale=0.1,mu=4",Size = dim(size_ClusterB_perm)[2],X= DataSize,conf = T,legendC = F)


plot_of_power(Data=size_ClusterC_perm,yval=c(0,1),xlab = "Pop. Size",title = "Clust: scale=0.5,mu=4",Size = dim(size_ClusterC_perm)[2],X=DataSize,conf = T,legendC = F)

plot_of_power(Data=size_pois_perm,yval=c(0,1),xlab = "Pop. Size",title = "Poisson",Size = dim(size_pois_perm)[2],X= DataSize,conf = T,legendC = F)

par(mfrow = c(1,1))

# We load the result for the LOF tests used on different pop. sizes. 
size_MaternA_OFT= readRDS(file = "Size_MaternA_OFT.Rdata") <=  0.05
size_MaternB_OFT= readRDS(file = "Size_MaternB_OFT.Rdata") <=  0.05
size_ClusterA_OFT= readRDS(file = "Size_ClusterA_OFT.Rdata") <=  0.05
size_ClusterB_OFT= readRDS(file = "Size_ClusterB_OFT.Rdata") <=  0.05
size_ClusterC_OFT= readRDS(file = "Size_ClusterC_OFT.Rdata") <=  0.05
size_pois_OFT= readRDS(file = "Size_pois_OFT.Rdata") <=  0.05


size_MaternA_OFNP= readRDS(file = "Size_MaternA_OFNP.Rdata") <=  0.05
size_MaternB_OFNP= readRDS(file = "Size_MaternB_OFNP.Rdata") <=  0.05
size_ClusterA_OFNP= readRDS(file = "Size_ClusterA_OFNP.Rdata") <=  0.05
size_ClusterB_OFNP= readRDS(file = "Size_ClusterB_OFNP.Rdata") <=  0.05
size_ClusterC_OFNP= readRDS(file = "Size_ClusterC_OFNP.Rdata") <=  0.05
size_pois_OFNP= readRDS(file = "Size_pois_OFNP.Rdata")  <=  0.05

size_MaternA_OFST= readRDS(file = "Size_MaternA_OFST.Rdata") <=  0.05
size_MaternB_OFST= readRDS(file = "Size_MaternB_OFST.Rdata") <=  0.05
size_ClusterA_OFST= readRDS(file = "Size_ClusterA_OFST.Rdata") <=  0.05
size_ClusterB_OFST= readRDS(file = "Size_ClusterB_OFST.Rdata") <=  0.05
size_ClusterC_OFST= readRDS(file = "Size_ClusterC_OFST.Rdata") <=  0.05
size_pois_OFST= readRDS(file = "Size_pois_OFST.Rdata") <=  0.05



# We make a function to turn convert the result into a latex table
table_of_LOFsize = function(MM, title= ""){
  MM = rbind(c("20","25","30","35","40" ) ,MM)
MM = data.frame(MM,row.names = c("Population Size","Hardcore II: r=0.05", "Hardcore II: r=0.02", "Cluster:scale=0.1,mu=1", "Cluster:scale=0.1,mu=4", "Cluster:scale=0.5,mu=4", "Poisson"))
tableMM = xtable(MM,auto=T)
print(tableMM,include.colnames = FALSE,
      hline.after = c(0,1,nrow(tableMM)))}

MM = rowmsd(size_MaternA_OFT)
MM = rbind(MM,rowmsd(size_MaternB_OFT))
MM = rbind(MM,rowmsd(size_ClusterA_OFT))
MM = rbind(MM,rowmsd(size_ClusterB_OFT))
MM = rbind(MM,rowmsd(size_ClusterC_OFT))
MM = rbind(MM,rowmsd(size_pois_OFT))

#MM = rbind(MM,c(msd(msd(size_MaternA_OFST)),msd(size_MaternB_OFST),msd(size_ClusterA_OFST),
#               msd(size_ClusterB_OFST),msd(size_ClusterC_OFST),msd(size_pois_OFST)))

#MM = rbind(MM,c(msd(size_MaternA_OFNP),msd(size_MaternB_OFNP),msd(size_ClusterA_OFNP),
#                msd(size_ClusterB_OFNP),msd(size_ClusterC_OFNP),msd(size_pois_OFNP)))


table_of_LOFsize(MM,title = "T-stat")

MM = rowmsd(size_MaternA_OFNP)
MM = rbind(MM,rowmsd(size_MaternB_OFNP))
MM = rbind(MM,rowmsd(size_ClusterA_OFNP))
MM = rbind(MM,rowmsd(size_ClusterB_OFNP))
MM = rbind(MM,rowmsd(size_ClusterC_OFNP))
MM = rbind(MM,rowmsd(size_pois_OFNP))
MM  = signif(MM,digits = 4)
table_of_LOFsize(MM,title = "Nearest Point")

MM = rowMeans(size_MaternA_OFST)
MM = rbind(MM,rowMeans(size_MaternB_OFST))
MM = rbind(MM,rowMeans(size_ClusterA_OFST))
MM = rbind(MM,rowMeans(size_ClusterB_OFST))
MM = rbind(MM,rowMeans(size_ClusterC_OFST))
MM = rbind(MM,rowMeans(size_pois_OFST))
MM  = signif(MM,digits = 4)
table_of_LOFsize(MM,title = "Spike-time")