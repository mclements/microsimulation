require(microsimulation)
library(snow)
library(Rmpi)

#mc.cores <- max(1, mpi.universe.size() - 1)
mc.cores <- mpi.universe.size()
cl <- makeMPIcluster(mc.cores)
cat(sprintf("Running with %d workers\n", length(cl)))
clusterEvalQ(cl, {library(microsimulation); NULL })

#print(test <- callFhcrc(1e6, mc.cores=1))
#print(test <- callFhcrc(1e6, mc.cores=7, cl=cl))
print(test <- callFhcrc(1e7, mc.cores=mc.cores, cl=cl))

stopCluster(cl)
mpi.quit()
