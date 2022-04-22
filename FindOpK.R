#We wish it fine the optimal K values for LOF when used on point patterns
library(spatstat)
library(utils)
library(ppMeasures)
#we define function there will be needed
stDistPP = function(X,Y,...){
  p1 = cbind(X$x,X$y)
  p2 = cbind(Y$x,Y$y)
  return(stDist(p1, p2, ...))
}

distMppp = function(X,nx=3,ny=nx,method=1,minpoints=20,sumfunc=Kest,...){
  
  n = length(X)
  M = matrix(0,n,n)
  q=1
  for (i in c(1:n)){
    for (j in c(q:n)){
      if (method==1 ){
        tempstore = OutlierPPP_Permu(X[[i]],X[[j]],nx=nx,ny=ny,minpoints=minpoints,use.tbar=1,sumfunc=sumfunc,...)
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
stDistPP(Data[[1]],Data[[2]],pm=1,bypassCheck=T)

R_10 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(10)
  Y = rpoint(10)
  Z = rpoint(10)
  Result =  stDistPP(X,Y,pm=1,bypassCheck=T)$distance +stDistPP(Y,Z,pm=1,bypassCheck=T)$distance - stDistPP (X,Z,pm=1,bypassCheck=T)$distance
  R_10[i] = Result < 0
  if (Result < 0){
    X10 = X
    Y10 = Y
    Z10 = Z
  }
}

R_3 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(3)
  Y = rpoint(3)
  Z = rpoint(3)
  Result =  stDistPP(X,Y,pm=1,bypassCheck=T)$distance +stDistPP(Y,Z,pm=1,bypassCheck=T)$distance - stDistPP (X,Z,pm=1,bypassCheck=T)$distance
  R_3[i] = Result < 0
  if (Result < 0){
    X3 = X
    Y3 = Y
    Z3 = Z
  }
}
sum(R_3)

R_100 = c(1:1000)
for (i in c(1:1000)){
  X = rpoint(100)
  Y = rpoint(100)
  Z = rpoint(100)
  Result =  stDistPP(X,Y,pm=1,bypassCheck=T)$distance +stDistPP(Y,Z,pm=1,bypassCheck=T)$distance - stDistPP (X,Z,pm=1,bypassCheck=T)$distance
  R_100[i] = Result < 0
  if (Result < 0){
    X3 = X
    Y3 = Y
    Z3 = Z
  }
}
sum(R_100)

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
  plot(kval,mean_val,main=title,ylim = c(min(min_val),max(max_val)),xlab="k",ylab="outlier factor LOF")
  lines(max_val, lty =2, col=2 )
  lines(min_val,lty =3,col=3)
  arrows(kval, mean_val-sdev, kval, mean_val+sdev, length=0.01, angle=90, code=3) #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
  }

# We make the plot with the spike distance function
OF_st = outlier_factors_PP(Data[1:100],method = 3,nx=2,ny=2,minpoints = 4,k = c(2:50),pm=1,bypassCheck=T)

plot_k(OF_st,c(2:50))
legend("right",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_st,2,sd)
plot(sdev)



# We make the plot with the nearest neighbor function 
OF_nn = outlier_factors_PP(Data[1:100],method = 2,nx=2,ny=2,minpoints = 4,k = c(2:50))

plot_k(OF_nn, c(2:50))
legend("topright",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))
sdev = apply(OF_nn,2,sd)
plot(sdev)

#https://stackoverflow.com/questions/10389967/common-legend-for-multiple-plots-in-r
par(oma = c(4,1,1,1), mfrow = c(1, 2), mar = c(1, 2, 2, 1))
plot_k(OF_st,c(2:50),title = "Spike-time")
#legend("right",cex =0.50, legend =c("max","mean with std","min"),lty=c(2,1,3))

plot_k(OF_nn, c(2:50),title = "Nearest neighbor")
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom",cex =0.75, legend =c("max","mean with std","min"),lty=c(2,1,3),col=c("red","black","green"),horiz = TRUE)

par(mfrow=c(1,1))




TestNorm = ppp(x=rnorm(100),y=rnorm(100),window = owin(c(-2,2),c(-2,2)))
M_Norm = crossdist(TestNorm,TestNorm)
k=c(2:50)
n = TestNorm$n
m = length(k)
Result = matrix(0,nrow = n,ncol = m)
colnames(Result) <- k
for (i in  c(1:m)){
  for (j in c(1:n)){
    Result[j,i] = l_outlier_factor(M_Norm,k[i],j)
  }}
plot_k(Result, c(2:50))

TestNorm2 = ppp(x=rnorm(400),y=rnorm(400),window = owin(c(-4,4),c(-4,4)))
M_Norm2 = crossdist(TestNorm2,TestNorm2)
k=c(2:50)
n = TestNorm2$n
m = length(k)
Result2 = matrix(0,nrow = n,ncol = m)
colnames(Result2) <- k
for (i in  c(1:m)){
  for (j in c(1:n)){
    Result2[j,i] = l_outlier_factor(M_Norm2,k[i],j)
  }}
plot_k(Result2, c(2:50))


Mst = distMppp(Data[1:100],method = 3,pm=1,bypassCheck=T)
loc1 <- cmdscale(Mst)
x1 <- loc1[, 1]
y1 <- loc1[, 2]
plot(x1, y1, xlab = "", ylab = "",main = "cmdscale(Mst)")

Mnn = distMppp(Data[1:100],method = 2)
loc2 <- cmdscale(Mnn)
x2 <- loc2[, 1]
y2 <- loc2[, 2] 

marks2 = OF_nn[,25]
marks2[marks2<1] = 1
ppnn = ppp(x=x2,y=y2,marks = log(marks2),window = owin(range(x2),range(y2)))
plot(ppnn,)
plot(x2, y2, xlab = "", ylab = "",main = "cmdscale(Mnn)")




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


#makes plot to show problem with small k value

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

