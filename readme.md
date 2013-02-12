Microsimulation package for R
=============================

This is a relatively simple microsimulation package for R.

The package includes an R implementation of discrete event simulation, building on several R5 classes. The package also includes a more extensive C++ library, building on Rcpp and several other C/C++ libraries.

At its core, the C++ implementation combines the
[SSIM](http://www.inf.usi.ch/carzaniga/ssim/index.html) C++ discrete event simulation library with the RngStream C library for random numbers. Rcpp is used to ease the burden of passing parameters between R and C/C++.

