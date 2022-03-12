#' Studentized Permutation Test for Two Samples of Functional Data
#' 
#' Perform a studentized permutation test of equal distribution (mean) for two 
#' samples of functional data.
#' @param foos1,foos2 arrays of real-valued data, of dimension \eqn{n x m_1} and 
#'   \eqn{n x m_2}, where \eqn{m_1}, \eqn{m_2} are the group size, and \eqn{n} 
#'   is the length of the function vectors.
#' @param use.tbar logical, defaults to \code{FALSE}. Whether to apply 
#'   studentization after integration, see `Details'.
#' @param nperm \code{NULL} or an integer giving the number of random 
#'   permutations. If \code{NULL}, all permutations are used. Only feasible if 
#'   group sizes \eqn{m_1}, \eqn{m_2} do not exceed 10. Default value is 25000,
#'   thus all permutations will be used on groups with \code{m_1 = m_2 =9}.
#' @details Test statistics are integral means of studentized square distances 
#'   between the group means. The test statistics are closely related, but not 
#'   equal to Hotelling's two-sample T-squared statistic. It is assumed that 
#'   the functions all have the same, equidistant arguments. Depending on the 
#'   value of \code{use.tbar}, the test statistic is either
#'   \itemize{
#'   \item \eqn{T = mean [ (mu_1(x)-mu_2(x))^2 / (s_1^2(x)/m_1 + s_2^2(x)/m_2) ]}    or
#'   \item \eqn{Tbar = mean [ (mu_1(x)-mu_2(x))^2 ] / mean[ s_1^2(x)/m_1 + s_2^2(x)/m_2 ]}
#'   }
#'   where \eqn{mu_1(x), mu_2(x)} and \eqn{s_1^2(x), s_2^2(x)} are within group 
#'   means and variances at a fixed argument \eqn{x}, and the overall
#'   the mean is taken over all \eqn{n} arguments \eqn{x}.
#'   
#'   If \code{nperm == NULL}, the exact test with all permutations is used
#'   (combinations, for symmetry reasons). This may cause memory or computing 
#'   time issues.
#'   If \code{nperm} is given as an integer, the permutations are \code{sample}d randomly, 
#'   unless \code{nperm} is larger than the number of disjoint combinations. In that
#'   case, the exact version is used. 
#' @return 
#' A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the test statistic,}
#' \item{p.value}{the p-value of the test,}
#' \item{alternative}{a character string describing the alternative hypothesis,}
#' \item{method}{a character string indicating what type of test was performed,}
#' \item{data.name}{a character string giving the name(s) of the data.}
#'    
#' @author Ute Hahn,  \email{ute@@imf.au.dk}
#' @export
#' @source Hahn(2012), with slight modification (using the mean instead of the 
#'    integral, in order to avoid having to pass the arguments of the functions)
#' @references Hahn, U. (2012) A Studentized Permutation Test for the Comparison 
#' of Spatial Point Patterns. \emph{Journal of the American Statistical 
#' Association},  \bold{107} (498), 754--764.
#' @keywords htest
#' @keywords robust
#' @keywords nonparametric
#' @keywords ts    

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



OutlierPPP_Permu = function(Outlier,PPP,nx,ny=nx,minpoints=20,use.tbar=FALSE,rinterval=NULL,nperm=999,sumfunc=Kest,...){
  grid1 = quadrats(Outlier,nx=nx,ny=ny)
  splitOutlier = split(Outlier,f=grid1)
  OutlierStat = NULL
  for(i in c(1:(nx*ny))){
    if(splitOutlier[[i]]$n >= minpoints){
      if (is.null(OutlierStat)){
        TEMPF =  sumfunc(splitOutlier[[i]],r=rinterval)
        if (is.null(rinterval)){
          rinterval = TEMPF$r
        }
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
rpoint(100,nsim=20)
OutlierPPP_Permu(Outlier = rpoint(100),PPP= rpoint(100,nsim=20),minpoints = 10,nx=3)

library(spatstat)

library(tictoc)

Outliertest1V1 = c(1:100)
Outliertest1V1result = c(1:100)
Outliertest1V2 = c(1:100)
Outliertest1V2result = c(1:100)
for (i in c(1:100)){
  Outlier = rpoint(100)
  PPP= rpoint(100,nsim=20)
  tic()
  Outliertest1V1result[i] = OutlierPPP(Outlier = Outlier,PPP= PPP,minpoints = 10,nx=3)$p.value
  T1 = toc()
  Outliertest1V1[i] =  T1$toc-T1$tic
  tic()
  Outliertest1V2result[i] = OutlierPPP_Permu(Outlier = Outlier,PPP= PPP,minpoints = 10,nx=3)$p.value
  T2 = toc()
  Outliertest1V2[i] =  T2$toc-T2$tic
  
}

mean(Outliertest1V1)
mean(Outliertest1V2)

Outliertest2V1 = c(1:100)
Outliertest2V1result = c(1:100)
Outliertest2V2= c(1:100)
Outliertest2V2result = c(1:100)
for (i in c(1:100)){
  Outlier = rMatClust(100,scale=0.1,mu=4)
  PPP= rpoint(100,nsim=20)
  tic()
  Outliertest2V1result[i] = OutlierPPP(Outlier = Outlier,PPP= PPP,minpoints = 10,nx=3)$p.value
  T1 = toc()
  Outliertest2V1[i] =  T1$toc-T1$tic
  tic()
  Outliertest2V2result[i] = OutlierPPP_Permu(Outlier = Outlier,PPP= PPP,minpoints = 10,nx=3)$p.value
  T2 = toc()
  Outliertest2V2[i] =  T2$toc-T2$tic
  
}

mean(Outliertest2V1)
mean(Outliertest2V2)
