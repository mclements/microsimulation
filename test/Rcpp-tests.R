require(Rcpp)
require(inline)
require(R.utils)

###########################
## Test Table<Key,Value> ##
###########################

testTable2 <- function(...) R.utils::intToBin(testTable(...))

src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
// [[Rcpp::export]]
SEXP testTable(DataFrame d, double x) {
  Table<double,double> lookup(d,"x","y");
  return wrap(lookup(x));
}
'
sourceCpp(code=src)
lookup <- data.frame(x=0:1,y=0:1)
testTable2(lookup,-1)
testTable2(lookup,0)
testTable2(lookup,0.1)
testTable2(lookup,2)

src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
using namespace std;
// [[Rcpp::export]]
SEXP testTable(DataFrame d, double x, double y) {
  Table<double,double,double> lookup(d,"x","y","z");
  return wrap(lookup(x,y));
}
'
sourceCpp(code=src)
lookup <- transform(expand.grid(data.frame(x=0:1,y=0:1)),z=0:3)
testTable2(lookup,-1,-1)
testTable2(lookup,0,0)
testTable2(lookup,1,0)
testTable2(lookup,0,1)
testTable2(lookup,2,2)

src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
#include <boost/tuple/tuple.hpp>
typedef boost::tuple<double,double,double> Tuple;
// [[Rcpp::export]]
SEXP testTable(DataFrame d, double x, double y, double z) {
  Table<double,double,double,double> lookup(d,"x","y","z","val");
  return wrap(lookup(x,y,z));
}
'
sourceCpp(code=src)
lookup <- transform(expand.grid(data.frame(x=0:1,y=0:1,z=0:1)),val=0:7)
testTable2(lookup,-1,-1,-1)
testTable2(lookup,0,0,0)
testTable2(lookup,1,0,1)
testTable2(lookup,0,1,2)
testTable2(lookup,2,2,2)

src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
#include <boost/tuple/tuple.hpp>
typedef boost::tuple<double,double,double,double> Tuple;
// [[Rcpp::export]]
SEXP testTable(DataFrame d, double x, double y, double z, double a) {
  Table<double,double,double,double,double> lookup(d,"x","y","z","a","val");
  return wrap(lookup(x,y,z,a));
}
'
sourceCpp(code=src)
lookup <- data.frame(expand.grid(data.frame(x=0:1,y=0:1,z=0:1,a=0:1)),val=0:15)
testTable2(lookup,-1,-1,-1,-1)
testTable2(lookup,0,0,0,0)
testTable2(lookup,1,0,1,1)
testTable2(lookup,0,1,2,2)
testTable2(lookup,2,2,2,3)

src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
#include <boost/tuple/tuple.hpp>
typedef boost::tuple<double,double,double,double,double> Tuple;
// [[Rcpp::export]]
SEXP testTable(DataFrame d, double x, double y, double z, double a, double b) {
  Table<double,double,double,double,double,double> lookup(d,"x","y","z","a","b","val");
  return wrap(lookup(x,y,z,a,b));
}
'
sourceCpp(code=src)
lookup <- data.frame(expand.grid(data.frame(x=0:1,y=0:1,z=0:1,a=0:1,b=0:1)),val=0:31)
testTable2(lookup,-1,-1,-1,-1,-1)
testTable2(lookup,0,0,0,0,0)
testTable2(lookup,1,0,1,1,0)
testTable2(lookup,0,1,2,2,2)
testTable2(lookup,2,2,2,3,3)


