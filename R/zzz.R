.onLoad <- function (lib, pkg) {
    library.dynam("microsimulation", pkg, lib)
    .microsimulation.init()
}

.onUnload <- function (libpath) {
  .microsimulation.exit()
  library.dynam.unload("microsimulation", libpath)
}
