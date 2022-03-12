
library(spatstat)

TESTFUNC = Kest
testpyramidal <- studpermu.test(pyramidal, Neurons ~ group ,summaryfunction = Fest)
testpyramidal$statistic


STU_TEST = studpermu.test(list(rpoint(100,nsim=1),rpoint(100,nsim=1)))
plot(STU_TEST)
STU_TEST$p.value


temp1 = rpoint(1000, nsim = 1)
temp2 = rpoint(1000, nsim = 1)

AA = quadrats(temp1,nx=3)

splitPP1 = split(temp1,f=AA)
splitPP2 = split(temp2,f=AA)


origon = function(X){
  DataWindow = X$window
  Data_xrange = DataWindow$xrange
  Data_yrange = DataWindow$yrange
  shift(X,c(-Data_xrange[1],-Data_yrange[1]))
}

newsplitPP1 = lapply(splitPP1,origon)
newsplitPP2 = lapply(splitPP2,origon)  


studpermu.test(list(splitPP1,splitPP2))
studpermu.test(list(newsplitPP1,newsplitPP2))

UteMethod = function(PPP1,PPP2,nx,ny=nx,rinterval=NULL,nperm=1,minpoints=10,sumfunc="Kest"){
  grid1 = quadrats(PPP1,nx=nx,ny=ny)
  grid2 = quadrats(PPP2,nx=nx,ny=ny)
  splitPP1 = split(PPP1,f=grid1)
  splitPP2 = split(PPP2,f=grid2)
  if (sumfunc == "Kest")
    return(studpermu.test(list(splitPP1,splitPP2),minpoints=minpoints,rinterval	
                          =rinterval,nperm=nperm, summaryfunction = Kest))
  else if (sumfunc =="Fest"){
    return(studpermu.test(list(splitPP1,splitPP2),minpoints=minpoints,rinterval	
                          =rinterval,nperm=nperm, summaryfunction = Fest))
  }
  else if (sumfunc =="Gest"){
    return(studpermu.test(list(splitPP1,splitPP2),minpoints=minpoints,rinterval	
                          =rinterval,nperm=nperm, summaryfunction = Gest))
  }
  else if (sumfunc =="Jest"){
    return(studpermu.test(list(splitPP1,splitPP2),minpoints=minpoints,rinterval	
                          =rinterval,nperm=nperm, summaryfunction = Jest))
  }
}

temp1 = rpoint(200, nsim = 2)
V = UteMethod(temp1[[1]],temp1[[2]],nx=3)
V$p.value
testL = rpoint(200, nsim = 4)

distMppp = function(X,nx,ny=ny,method=1,minpoints=20){
  pvaluetrans = function(p,method){
    if (method==1){
      return(1/p -1)
    }
    else if (method==2){
      return(1-p)
    }
  }
  n = length(X)
M = matrix(0,n,n)
q=1
for (i in c(1:n)){
  for (j in c(q:n)){
    tempstore=UteMethod(X[[i]],X[[j]],nx=3,minpoints=minpoints)
    
    M[i,j]=tempstore$statistic
    M[j,i]=M[i,j]
  }
  q = q+1
}
  return(M)
}
library(cluster)

M1 = distMppp(testL,nx=3,method = 1)
M2 = distMppp(testL,nx=3,method = 2)
hctest1= agnes(M1)
hctest2= agnes(M2)
plot(hctest1)
plot(hctest2)

testfortrekant = function(M){
  temp=list()
  n = dim(M)[1]
  for (i in c(1:n)){
    for (j in c(1:n)){
      for (l in c(1:n)){
        if (M[i,j] > M[i,l] + M[l,j]){
          temp = append(temp,list(c(i,j,l)))
        }
      }}}
  return(temp)
}
testfortrekant(M1)


Pois1 = rpoispp(200,nsim=10)
Matern1 = rMaternI(250,r=0.01,nsim=10)
Clust1 = rMatClust(40,scale=0.05,mu=5,nsim = 10)
L = c(Pois1,Matern1,Clust1)
MDtest1 = distMppp(L,nx=3,method = 1)

hctestPP1= agnes(MDtest1,method = "average")
plot(hctestPP1)

hctestPP1= agnes(MDtest1^2)
plot(hctestPP1)

hctestPP2 = diana(MDtest1)
plot(hctestPP2)
library(colorspace)
library(dendextend)

make_dend = function(hc,k,label){
  n = length( hc$order)
  dend <- as.dendrogram(hc)
  dend <- rotate(dend, 1:n)
  dend <- color_branches(dend, k=k)
  
  labels_colors(dend) <-
    rainbow_hcl(k)[sort_levels_values(
      as.numeric(as.factor(label))[order.dendrogram(dend)]
    )]
  
  labels(dend) <- paste(as.character(label)[order.dendrogram(dend)],
                        "(",labels(dend),")", 
                        sep = "")
  #dend <- hang.dendrogram(dend,hang_height=0.1)
  dend <- set(dend, "labels_cex", 0.5)
  return(dend)
}
label = c(rep("Pois",10),rep("Matern1",10),rep("Clust",10))
dend = make_dend(hctestPP1,k=3,label)

plot(dend,horiz =  TRUE,  nodePar = list(cex = .007))

Pois2 = rpoispp(200,nsim=5)
Matern2 = rMaternI(250,r=0.03,nsim=5)
MDtest2 = distMppp(c(Pois2,Matern2),nx=2,method = 1,minpoints=10)
hctestPP3= agnes(MDtest2^2,method = "average")
plot(hctestPP3)

windowEX= owin(xrange = c(0,3),yrange=c(0,3))
Pois3 = rpoispp(200,nsim=5,win=windowEX)
Matern3 = rMaternI(250,r=0.03,nsim=5,win=windowEX)
Clust3 = rMatClust(40,scale=0.05,mu=5,nsim = 5,win=windowEX)
L = c(Pois1,Matern1,Clust1)
MDtest3 = distMppp(L,nx=5,method = 1,minpoints=10)

hctestPP3= agnes(MDtest3,method = "average")
plot(hctestPP3)



OutlierPPP = function(Outlier,PPP,nx,ny=nx, comOregon=0,minpoints=20,sumfunc = "Kest",rinterval=NULL,nperm=999){
  grid1 = quadrats(Outlier,nx=nx,ny=ny)
  splitOutlier = split(Outlier,f=grid1)
  n = length(PPP)
  splitPP= list()
  for (i in c(1:n)){
    grid2 = quadrats(PPP[[i]],nx=nx,ny=ny)
    splitPP = append(splitPP,split(PPP[[i]],f=grid2))
  }

  if (comOregon==0){
    if (sumfunc == "Kest")
    return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
=rinterval,nperm=nperm, summaryfunction = Kest))
     else if (sumfunc =="Fest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Fest))
     }
    else if (sumfunc =="Gest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Gest))
    }
    else if (sumfunc =="Jest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Jest))
    }
  } else{
    newsplitPP1 = lapply(splitOutlier,origon)
    newsplitPP2 = lapply(splitPP,origon)  
    if (sumfunc == "Kest")
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Kest))
    else if (sumfunc =="Fest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Fest))
    }
    else if (sumfunc =="Gest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Gest))
    }
    else if (sumfunc =="Jest"){
      return(studpermu.test(list(splitOutlier,splitPP),minpoints=minpoints,rinterval	
                            =rinterval,nperm=nperm, summaryfunction = Jest))
    }
  }
  
}


pois4 = rpoispp(100,nsim=20)
Matern4a = rMaternI(100,r=0.05,nsim=1)
Matern4b = rMaternI(100,r=0.02,nsim=1)
Clust4a = rMatClust(100,scale=0.1,mu=1)
Clust4b = rMatClust(100,scale=0.1,mu=4)
Clust4c = rMatClust(100,scale=0.05,mu=4)


OutlierPPP(Outlier = rpoispp(100), PPP = pois4, nx=2,minpoints = 10)

QQQ = OutlierPPP(Outlier = rpoispp(100), PPP = pois4, nx=2,minpoints = 10)


Resultmatrix = matrix(0,nrow=2,ncol = 5)




OutlierPPP(Outlier = Matern4a, PPP = pois4, nx=2,minpoints = 5)$p.value
OutlierPPP(Outlier = Matern4b, PPP = pois4, nx=2,minpoints = 5)$p.value
OutlierPPP(Outlier = Clust4a, PPP = pois4, nx=2,minpoints = 5)$p.value
OutlierPPP(Outlier = Clust4b, PPP = pois4, nx=2,minpoints = 5)$p.value
OutlierPPP(Outlier = Clust4c, PPP = pois4, nx=2,minpoints = 5)$p.value

sumfuncL = c("Kest","Fest","Gest")
ppL = list(Matern4a=Matern4a,Matern4b=Matern4b,Clust4a=Clust4a,Clust4b=Clust4b,Clust4c=Clust4c)
plot(as.ppplist(ppL))
Resultmatrix = matrix(0,nrow = length(sumfuncL), ncol = length(ppL))
rownames(Resultmatrix) <- sumfuncL
colnames(Resultmatrix) <- names(ppL)
for (i in c(1:length(sumfuncL))) {
  for(j in c(1:length(ppL)))
    Resultmatrix[i,j] = OutlierPPP(Outlier = ppL[[j]], PPP = pois4, nx=2,minpoints = 5, sumfunc = sumfuncL[i],r=c(0,0.000977))$p.value
}
Resultmatrix

OutlierPPP(Outlier = Matern4b, PPP = pois4, nx=2,minpoints = 5, sumfunc = "Fest")






