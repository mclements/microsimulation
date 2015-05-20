require(microsimulation)
for (i in c(1, 2, 4, 8)){
  cat(paste0("Number of cores:", i,"\n"))
  callFhcrc(n = 1e7, mc.cores = i)
}

