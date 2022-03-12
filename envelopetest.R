library(spatstat)
library(GET)
X <- spruces

pois4env  = rpoispp(100,nsim=20)
Matern4aenv = rMaternI(100,r=0.05,nsim=1)
Matern4benv  = rMaternI(100,r=0.02,nsim=1)
Clust4aenv  = rMatClust(100,scale=0.1,mu=1)
Clust4benv  = rMatClust(100,scale=0.1,mu=4)
Clust4cenv  = rMatClust(100,scale=0.05,mu=4)


env <-  envelope(Matern4aenv, nsim=999, savefuns=TRUE,
                  simulate=expression(runifpoint(ex=X)),
                  verbose=FALSE)
envelope(rpoispp(100,nsim=10), nsim=999, savefuns=TRUE,
         simulate=expression(runifpoint(ex=X)),
         verbose=FALSE)
env2 <- envelope(Matern4aenv, nsim=length(pois4env), savefuns=TRUE,
         simulate=pois4env,
         verbose=FALSE)
test1_erl = global_envelope_test(env2,type = "erl")
attributes(test1_erl)$p
global_envelope_test(env,type = "erl")
test1_rank = global_envelope_test(env2,type = "rank")

env2 <- envelope(Matern4aenv, nsim=length(pois4env), savefuns=TRUE,
                 simulate=pois4env,
                 verbose=FALSE)

frank.fanova(env,nsim = 100)

testenvelope = GET.necdf(env,nsim=100)
plot(testenvelope)
