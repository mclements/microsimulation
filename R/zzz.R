
#' Internal function
#' @param lib library string
#' @param pkg package string
#' @rdname Utilities
.onLoad <- function (lib, pkg) {
#  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
#   invisible(rnorm(1)) # kludge to initialise the random seed
  microsimulation.init()
}

#' Internal function
#' @param libpath library path string
#' @rdname Utilities
.onUnload <- function (libpath) {
  microsimulation.exit()
  library.dynam.unload("microsimulation", libpath)
}
