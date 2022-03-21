library(spatstat)

TP = ppp(x=c(0.2,0.2,0.2,0.4,0.3),y=c(0.8,0.9,0.7,0.8,0.2))
plot(TP)
text(TP , labels= c("x","x","x","z","y"),pos=4)

X= matrix(c(0.2,0.8,0.2,0.9,0.2,0.7),byrow = T,3,2)

dXY = sum(sqrt(rowSums((X-c(0.3,0.2))^2)))
dYX = sum(sqrt(((X[3,]-c(0.3,0.2))^2)))

dXZ = sum(sqrt(rowSums((X-c(0.4,0.8))^2)))
dZX = sum(sqrt(((X[2,]-c(0.4,0.8))^2)))


dZY = sum(sqrt((c(0.3,0.2)-c(0.4,0.8))^2))
dYZ = sum(sqrt((c(0.3,0.2)-c(0.4,0.8))^2))
  dXZ+ dZX + dZY + dYZ -(dXY+ dYX)

  X2= matrix(c(0.2,0.8,0.2,0.9,0.2,0.7 ,0.1,0.7 ),byrow = T,4,2)

  
  dXY = sum(sqrt(rowSums((X2-c(0.3,0.2))^2)))
  dYX = sum(sqrt(((X2[3,]-c(0.3,0.2))^2)))
  
  dXZ = sum(sqrt(rowSums((X2-c(0.4,0.8))^2)))
  dZX = sum(sqrt(((X2[2,]-c(0.4,0.8))^2)))
  
  
  dZY = sum(sqrt((c(0.3,0.2)-c(0.4,0.8))^2))
  dYZ = sum(sqrt((c(0.3,0.2)-c(0.4,0.8))^2))
  dXZ+ dZX + dZY + dYZ >(dXY+ dYX)
  
  
  X3= matrix( c(XR$x,XR$y),byrow = F,3,2)
  YRc = c(YR$x,YR$y)
  ZRc =  c(ZR$x,ZR$y)
  
  dXY = sum(sqrt(rowSums((X3-YRc)^2)))
  dYX = sqrt(min(rowSums((X3-YRc)^2)))
  
  dXZ =  sum(sqrt(rowSums((X3-ZRc)^2)))
  dZX = sqrt(min(rowSums((X3-ZRc)^2)))
  
  
  dZY = sum(sqrt((YRc-ZRc)^2))
  dYZ = sum(sqrt((YRc-ZRc)^2))
  dXZ+ dZX + dZY + dYZ > (dXY+ dYX)
  
  
  
  LX = coords(XR)
  LY = coords(YR)
  FF = function(x){
    sqrt(x[1]^2 +x[2]^2)
  }
  result = 0
  n = XR$n
  for (i in c(1:n)){
    result = result+  min(apply(LY - c(LX[i,1:2]),1 , FUN = FF))
  }
  return(result)
  
  
  LY - c(LX[2,1:2])
  (X3-YRc)[2,]
  LY - c(X3[2,1:2])
