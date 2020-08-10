

.microsimulationLdFlags <- function(){
    paste( '-L"', system.file( "lib", package = "microsimulation" ), '" -lmicrosimulation',
          sep = "" )
}

inlineCxxPlugin <- function(...) {
    ismacos <- Sys.info()[["sysname"]] == "Darwin"
    openmpflag <- if (ismacos) "" else "$(SHLIB_OPENMP_CFLAGS)"
    plugin <- Rcpp::Rcpp.plugin.maker(include.before = "#include <microsimulation.h>",
                                      libs           = paste(openmpflag,
                                                             .microsimulationLdFlags(),
                                                             "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),
                                      package        = "microsimulation")
    settings <- plugin()
    settings$env$PKG_CPPFLAGS <- paste("-I../inst/include", openmpflag)
    ## if (!ismacos) settings$env$USE_CXX11 <- "yes"
    settings
}
