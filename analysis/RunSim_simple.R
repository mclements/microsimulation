require(microsimulation)
for (i in c(1, 2, 4, 8)){
  print(paste("Number of cores:", i))
  callFhcrc(n = 1e7, mc.cores = i)
}

