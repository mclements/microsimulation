.onLoad <- function (lib, pkg) {
    library.dynam("rlecuyer", pkg, lib)
    .lec.init()
}

.onUnload <- function (libpath) {
  .lec.exit()
  library.dynam.unload("rlecuyer", libpath)
}
