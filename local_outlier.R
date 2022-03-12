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
  Result = matrix(0,nrow = n,ncol = m)
  colnames(Result) <- k
  for (i in  c(1:m)){
  for (j in c(1:n)){
    Result[j,i] = l_outlier_factor(M,k[i],j)
  }}
  return(Result)
}
pois_of  = rpoispp(100,nsim=20)
Matern_a_of = rMaternI(100,r=0.05,nsim=1)
Matern_b_of  = rMaternI(100,r=0.02,nsim=1)
Clust_a_of  = rMatClust(100,scale=0.1,mu=1)
Clust_b_of  = rMatClust(100,scale=0.1,mu=4)
Clust_c_of = rMatClust(100,scale=0.05,mu=4)
TEMP_PP =append(pois_of , c(list(Matern_a_of),list(Matern_b_of),list(Clust_a_of ),list(Clust_b_of ),list(Clust_c_of )))
outlier_factors_PP(TEMP_PP,k=3,nx=3,minpoints = 5)

outlier_factors_PP(TEMP_PP,k=6,nx=3,minpoints = 5)

TEMP_PP =append(rpoispp(100,nsim=40) , c(list(Matern_a_of),list(Matern_b_of),list(Clust_a_of ),list(Clust_b_of ),list(Clust_c_of )))
outlier_factors_PP(TEMP_PP,k=3,nx=3,minpoints = 5)

apply(outlier_factors_PP(TEMP_PP,k=c(10:20),nx=3,minpoints = 5),MARGIN = 1, FUN = max)
