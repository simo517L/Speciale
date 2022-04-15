# Simplified implementation of Simons experiment
# Ute, 13.4.22
#
# simulation experiment setup:
# check one candidate for outlier against a population of size 20
# by a studentized permutation test, using estimates of the K-function
# Simon has stored the candidates and the populations in a data file


#------------ preliminaries and setup ------------

library(spatstat)
source("permutationtest_pval.R")

# summary function: K-function with argument r = rvalues
# the choice of r-values follows an old recommendation of
# Ripley (?), to look at an interval from 0 to 1.25/sqrt(lambda)
# Here, we had lambda = 100
# we use isotropic edge correktion
# We want to discard patterns with too few points, less than minpoints.
# here it is done by marking the results as NA. Then we discard them later

sumfun_K <- function(pp, 
                     minpoints = 10, 
                     rvalues = seq(0, .125, 0.005)){
  if (npoints(pp) < minpoints) rep(NA, length(rvalues)) 
  else Kest(pp, r = rvalues, correction = "iso")$iso
}

# the experiments requires estimating a summary function on subwindows
# of the patterns (quadrats)
# Setup for the quadrats, minpoints, and r-values 
# (although already included as defaults in sumfun_K)

setup <- list(
  quadsizes = list(
    "1x3" = list(nx = 1, ny = 3),
    "2x2" = list(nx = 2, ny = 2),
    "2x3" = list(nx = 2, ny = 3),
    "3x3" = list(nx = 3, ny = 3),
    "3x4" = list(nx = 3, ny = 4),
    "4x4" = list(nx = 4, ny = 4)),
  minpoints = 10,
  rvalues = seq(0, .125, 0.005),
  sumfun = sumfun_K,
  nperm = 1000 # number of permutations in studpermute.test
  # chose a small value for pilot experiments. 
  # Could be larger when run on server
)


#---------- helper functions -----------------

# apply a summary function to a  point pattern,
# splitted into subquadrats
# returns a matrix

sumfun_on_quads <- function(pp, 
                            nx = 2, ny = nx,
                            minpoints = 10,
                            rvalues = seq(0, .125, 0.005),
                            sumfun = sumfun_K)
{
  quads <- quadrats(pp, nx, ny)
  pps <- split(pp, quads)
  result <- sapply(pps, sumfun, 
                   minpoints = minpoints,
                   rvalues = rvalues)
}

# apply sumfun_on_quads on all patterns in a list of point patterns.
# then remove results from quadrats that did not have enough points.

sumfun_array <- function(list_of_patterns,
                         nx = 2, ny = nx,
                         minpoints = 10,
                         rvalues = seq(0, .125, 0.005),
                         sumfun = sumfun_K)
{ # get a list of results for each pattern
  sumfunlist <- lapply(list_of_patterns, sumfun_on_quads,
                   nx = nx, ny = ny, 
                   minpoints = minpoints, rvalues = rvalues,
                   sumfun = sumfun)
  # transform this list into a matrix 
  # the individual results are matrices of dimension nr x nquad
  # we want to append all these matrices into one of dimension nr x (npp*nquad)
  npp <-   length(list_of_patterns) # number of patterns
  nr <-    length(rvalues) 
  nquad <- nx * ny # number of quadrats per pattern
  
  result <- array(0, c(nr, npp*nquad))
  for (i in 1:npp) 
    result[, (i-1)*nquad + (1:nquad)] <- sumfunlist[[i]]
  result
}



#----------- main program ---------

populations <- readRDS("DataPPP.Rdata")
candidates <- readRDS("poisPPP.Rdata")

ncandidates <- length(candidates)
npopulations <- length(populations)

popsize <- npopulations %/% ncandidates # should be 20

# get pvalues for candidates in i_candidates = 1: ncandidates
# For debugging purposes, this can be changed to a smaller number

icandidates <- 1 : ncandidates
nquadsizes = length(setup$quadsizes)

pvalues <- array(NA, c(nquadsizes, length(icandidates)))
# record also the number of quadrats that had enough points
nquads_cand <- nquads_pop <- array(0, c(nquadsizes, length(icandidates)))
row.names(pvalues) <- row.names(nquads_cand) <- row.names(nquads_pop) <- names(setup$quadsizes)

for (i in seq_along(icandidates)) {
  # get p-values for one candidate and the population it is tested against
  icand <- icandidates[i]
  candidate <- candidates[[icand]]
  population <- populations[popsize * (icand-1) + (1 : popsize)]
  
  # go through all quadrat sizes
  for (j in seq_along(setup$quadsizes)){
    qs <- setup$quadsizes[[j]]
    foos1 <- sumfun_on_quads(candidate, 
                             nx = qs$nx, ny = qs$ny,
                             minpoints = setup$minpoints,
                             rvalues = setup$rvalues,
                             sumfun = setup$sumfun)
    # keep only those that have enough points
    keep1 <- !is.na(foos1[1,])
    nq1 <- nquads_cand[j, i] <- sum(keep1)
    if (nq1 > 1) {
      foos2 <- sumfun_array(population,
                          nx = qs$nx, ny = qs$ny,
                          minpoints = setup$minpoints,
                          rvalues = setup$rvalues,
                          sumfun = setup$sumfun)
      keep2 <- !is.na(foos2[1,])
      nq2 <- nquads_pop[j, i] <- sum(keep2)
      if (nq2 > 1)
      pvalues[j, i] <- studpermute.pval(foos1[, keep1],
                                      foos2[, keep2], 
                                      nperm = setup$nperm)
    }  
  }
}

save(pvalues, nquads_cand, nquads_pop, 
     file = "Results_Permutationtest_Pois.RData")

##------------ checking uniformity ------------
# a qqplot
# plot(ppoints(length(icandidates)), sort(pvalues["2x2",]))
# abline(0,1)
