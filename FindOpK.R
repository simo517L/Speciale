#We wish it fine the optimal K values for LOF when used on point patterns


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

plot_k = function(X,kval){
  max_val = do.call(pmax,as.data.frame(t(X)))
  min_val = do.call(pmin,as.data.frame(t(X)))
  mean_val = colMeans(X)
  sdev = apply(X,2,sd)
  plot(kval,mean_val,ylim = c(min(min_val),max(max_val)))
  lines(max_val)
  lines(min_val)
  arrows(kval, mean_val-sdev, kval, mean_val+sdev, length=0.01, angle=90, code=3) #https://stackoverflow.com/questions/13032777/scatter-plot-with-error-bars
}


OF_st = outlier_factors_PP(Data[1:100],method = 3,nx=2,ny=2,minpoints = 4,k = c(1:50),pm=1,bypassCheck=T)


plot_k(OF_st[,-1],c(2:50))
sdev = apply(OF_st[,-1],2,sd)
plot(sdev)



OF_nn = outlier_factors_PP(Data[1:100],method = 2,nx=2,ny=2,minpoints = 4,k = c(2:50))

plot_k(OF_nn, c(2:50))

sdev = apply(OF_nn,2,sd)
plot(sdev)

TestNorm = ppp(x=rnorm(300),y=rnorm(300),window = owin(c(-3,3),c(-3,3)))
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
plot_k(Result[,1:14], c(1:14))
