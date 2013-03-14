.onLoad <- function (lib, pkg) {
  ## if (!exists(".Random.seed", envir = .GlobalEnv, 
  ##             inherits = FALSE))
  ##   invisible(rnorm(1)) # kludge to initialise the random seed
  library.dynam("microsimulation", pkg, lib)
  .microsimulation.init()
}

.onUnload <- function (libpath) {
  .microsimulation.exit()
  library.dynam.unload("microsimulation", libpath)
}
