# corrected version of permutation test, only p=value
# 14.4.22 uh
# Simon has a more sophisticated version

studpermute.pval <- function (foos1, foos2, use.tbar=FALSE, nperm = 1000)
{
  ##### preparations ----------------
  n <- dim(foos1)[1]
  m1 <- dim(foos1)[2]
  m2 <- dim(foos2)[2]
  
  # need at least 2 per group, and the dimensions have to be compatible
  if (m1 < 3 | m2 < 3 | dim(foos2)[1] != n){
    return(NA)
  }
  
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
    ss1 <- (SXX1 - (SX1^2 / m1)) / ((m1-1) * m1)
    ss2 <- (SXX2 - (SX2^2 / m2)) / ((m2-1) * m2)
    
    if (use.tbar) return (sum((mu1 -mu2)^2) / sum((ss1 + ss2))) else 
      return (mean((mu1 -mu2)^2 / (ss1 + ss2), na.rm=T))
  }
  
  Tvals <- apply(index1, 2, Tstatistic)
  # p-value  
  mean(Tvals >= Tvals[1])           
}

