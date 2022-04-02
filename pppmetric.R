library(spatstat)
library(dplyr)
X = rpoint(10)
Y = rpoint(10)

LX = coords(X)
LY = coords(Y)
LY[1:5,] = LX[1:5,]
LX
LY

setdiff(LY,LX)
prod(L[1,1:2] == L)

spike_timedist = function(X,Y,pd=1,pa=pd,pm=1){

  nx = X$n
  ny= Y$n
  LX = coords(X)
  LY = coords(Y)
  FF = function(x){
    sqrt(x[1]^2 +x[2]^2)
  }
  

  equalM = function(L,M){
    return(dim(M) == dim(L) && all(M == L))
  }
  
  
  tempF= function(M,result=0,L=list(),RV = c()){
    n = length(L)
    for (i in c(1:n)){
      if (n > 0 &&dim(M) == dim(L[[i]]) && all(M == L[[i]])){
        return(list(RV[i] + result,L,RV))
      } 
    }
    if (dim(M)[1] ==1 & dim(M)[2] ==1){
      result = result+M[1]*pm
      return(list(result,append(L,list(M)),c(RV,result)))
    } else if (dim(M)[1] ==1 & dim(M)[2] >1){
      result = result+min(M[1,])*pm + (dim(M)[2]-1)*pa 
      return(list(result,append(L,list(M)),c(RV,result)))
    }else if (dim(M)[1] >1 & dim(M)[2] == 1){
      result = result+min(M[,1])*pm + (dim(M)[1]-1)*pd
      return(list(result,append(L,list(M)),c(RV,result)))
    } else{
      temp = c(1:dim(M)[1])
      for (i in c(1:dim(M)[1])){
        X = min(M[i,])
        tempResult = result + X*pm
        newM = M[-i,-i]
        if (is.matrix(newM)==T){
          tempRL =  tempF(M=newM , tempResult,L=L,RV=RV)
          
          temp[i] =  tempRL[[1]]
          if (length(L) == 0){
            L = append(L,list(newM))
            RV = c(RV,tempRL[[1]])
          }
          else if  ( sum(sapply(L, equalM , M = newM))==0){
            L = append(L,list(newM))
            RV = c(RV,tempRL[[1]])
          }
          
        } else{
          newM = matrix(newM, nrow = dim(M)[1]-1, ncol =dim(M)[2]-1)    
          tempRL =  tempF(M=newM , tempResult,L=L,RV=RV)
          
          temp[i] =  tempRL[[1]]
          if (length(L) == 0){
            L = append(L,list(newM))
            RV = c(RV,tempRL[[1]])
          }
          else if  ( sum(sapply(L, equalM , M = newM))==0){
            L = append(L,list(newM))
            RV = c(RV,tempRL[[1]])
          }
        }
        
      }
      return(list(min(temp),L,RV))
    }
  }
  
  distM = apply(LY - c(LX[1,1:2]),1 , FUN = FF)
  for (i in c(2:nx)){
   distM = rbind(distM, apply(LY - c(LX[i,1:2]),1 , FUN = FF))
  }

  return(tempF(distM)[[1]])
}

spike_timedist2 = function(X,Y,pd=1,pa=pd,pm=1,alg = c("default", "IMA", "VP97", "MSU"), euclid = NULL,
lossOrder = 1, maxBranch = 4, eps = 10^-10, bypassCheck = FALSE){
  p1 = cbind(X$x,X$y)
  p2 = cbind(Y$x,Y$y)
  hold = stDist(p1, p2,pd=pd,pa=pa,pm=pm,alg = alg, euclid = euclid,
                lossOrder = lossOrder, maxBranch = maxBranch, eps = eps, 
                bypassCheck = bypassCheck)
  return(hold)
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
X = rpoint(3)
Y = rpoint(3)
Z = rpoint(3)
spike_timedist(X,Y)
nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)

for (i in c(1:1000)){
  X = rpoint(3)
  Y = rpoint(1)
  Z = rpoint(1)
  Result =  nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)
  if(Result < 0){
    XR = X
    YR = Y
    ZR = Z
  }
}

XR = ppp(x=c(0.2,0.2,0.2,0.1,0.1,0.1),y=c(0.8,0.9,0.7,0.8,0.9,0.7))
YR = ppp(x=c(0.4),y=c(0.8))
ZR = ppp(x=c(0.3),y=c(0.2))
marks(XR) <- factor("x")
marks(YR) <- factor("y")
marks(ZR) <- factor("z")

plot(superimpose(XR,YR,ZR), main = "")

nearest_point_metric(XR,YR) + nearest_point_metric(YR,ZR) - nearest_point_metric(XR,ZR)



for (i in c(1:1000)){
  print(i)
  X = rpoint(100)
  Y = rpoint(100)
  Z = rpoint(100)
  MM = distMppp(list(X,Y,Z),nx=2,ny=2,method=1,minpoints=20,sumfunc=Kest)
  Result =  MM[1,2] + M[2,3] - M[1,3]
  if(Result < 0){
    XRL = X
    YRL = Y
    ZRL = Z
  }
}
MM = distMppp(list(XRL,YRL,ZRL),nx=2,ny=2,method=1,minpoints=20,sumfunc=Kest)
MM[1,2] + M[2,3] - M[1,3]
Counter_Example2 = hyperframe(list(XRL,YRL,ZRL),row.names=c("X","Y","Z"))
plot(Counter_Example2,main = "")
names(Counter_Example2[1])
R_10 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(10)
  Y = rpoint(10)
  Z = rpoint(10)
  Result =  nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)
  R_10[i] = Result < 0
}

R_100 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(100)
  Y = rpoint(100)
  Z = rpoint(100)
  Result =  nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)
  R_100[i] = Result < 0
}

R_1 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(1)
  Y = rpoint(1)
  Z = rpoint(1)
  Result =  nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)
  R_1[i] = Result < 0
}


R_2= c(1:1000)
for (i in c(1:1000)){
  X = rpoint(2)
  Y = rpoint(2)
  Z = rpoint(2)
  Result =  nearest_point_metric(X,Y) + nearest_point_metric(Y,Z) - nearest_point_metric(X,Z)
  R_2[i] = Result < 0
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
      } else if (method ==3){
        tempstore = spike_timedist2(X[[i]],X[[j]], bypassCheck =T)
        M[i,j] = tempstore$distance
      }
      
      M[j,i]=M[i,j]
    }
    q = q+1
  }
  return(M)
}
library(cluster)

pois5 = rpoispp(100,nsim=5)
Matern5a = rMaternI(100,r=0.05,nsim=5)
Matern5b = rMaternI(100,r=0.02,nsim=5)
Clust5a = rMatClust(100,scale=0.1,mu=1,nsim=5)
Clust5b = rMatClust(100,scale=0.1,mu=4,nsim=5)
Clust5c = rMatClust(100,scale=0.05,mu=4,nsim=5)

L = c(pois5,Matern5a,Clust5b)
MDtest3 = distMppp(L,nx=3,method = 2)
hctestPP3= agnes(MDtest3)
plot(hctestPP3)

L1 = c(pois5,Matern5b,Clust5a)
MDtest5 = distMppp(L1,nx=3,method = 2)
hctestPP5= agnes(MDtest5)
plot(hctestPP5)


#DONT RUN
MDtest4 = distMppp(L,nx=3,method = 3)
hctestPP4= agnes(MDtest4)
plot(hctestPP4)

MDtest5 = distMppp(L1,nx=3,method = 3)
hctestPP5= agnes(MDtest5)
plot(hctestPP5)

