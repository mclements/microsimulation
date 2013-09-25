#ifndef RCPP_TABLE_H
#define RCPP_TABLE_H

#include <Rcpp.h>

#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

using namespace std;
using namespace Rcpp;

class Interpolate {
public:
  vector<double> x, y, slope;
  Interpolate(vector<double> inx, vector<double> iny) : 
    x(inx), y(iny) { 
    // calculate the slope between points
    for (size_t i=0; i<x.size()-1; i++) {
      slope.push_back((y[i+1]-y[i]) / (x[i+1]-x[i]));
    };
  }
  Interpolate(DataFrame df) {
    NumericVector x = df(0);
    NumericVector y = df(1);
    // calculate the slope between points
    for (int i=0; i<x.size()-1; i++) {
      slope.push_back((y[i+1]-y[i]) / (x[i+1]-x[i]));
    };
  }
  double approx(double xfind) {
    int i;
    if (xfind<=x[0]) return y[0];
    else if (xfind>=*(--x.end())) return *(--y.end());
    else {
      i = lower_bound(x.begin(), x.end(), xfind) - x.begin();
      return y[i]+slope[i]*(xfind-x[i]);
    }
  }
  double operator()(double xfind) {
    if (xfind<=x[0]) return y[0];
    int i = lower_bound(x.begin(), x.end(), xfind) - x.begin();
    return y[--i];
  }
};

template<class Index, class Outcome>
class Table1D {
public:
  typedef map<Index,Outcome,greater<Index> > Map;
  Map data;
  Table1D(DataFrame df) { 
    Vector<Rcpp::traits::r_sexptype_traits<Index>::rtype> df0 = df(0);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df1 = df(1);
    for (int i=0; i<df0.size(); i++) {
      data[df0[i]] = df1[i];
    }
  }
  virtual Outcome lookup(Index index) {
    return data.lower_bound(index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
};

template<class I1, class I2, class Outcome>
class Table2D {
public:
  typedef pair<I1,I2> Index;
  typedef map<Index,Outcome,greater<Index> > Map;
  Map data;
  Table2D(DataFrame df) { 
    Vector<Rcpp::traits::r_sexptype_traits<I1>::rtype> df0 = df(0);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df1 = df(1);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df2 = df(2);
    for (int i=0; i<df0.size(); i++) {
      data[make_pair(df0[i],df1[i])] = df2[i];
    }
  }
  virtual Outcome lookup(Index index) {
    return data.lower_bound(index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2) {
    return data.lower_bound(make_pair(i1,i2))->second;
  }
  virtual Outcome operator()(I1 i1, I2 i2) {
    return lookup(i1,i2);
  }
};

template<class I1, class I2, class I3, class Outcome>
class Table3D {
public:
  typedef boost::tuple<I1,I2,I3> Index;
  typedef map<Index,Outcome,greater<Index> > Map;
  Map data;
  Table3D(DataFrame df) { 
    Vector<Rcpp::traits::r_sexptype_traits<I1>::rtype> df0 = df(0);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df1 = df(1);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df2 = df(2);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df3 = df(3);
    for (int i=0; i<df0.size(); i++) {
      data[Index(df0[i],df1[i],df2[i])] = df3[i];
    }
  }
  virtual Outcome lookup(Index index) {
    return data.lower_bound(index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2, I3 i3) {
    return data.lower_bound(Index(i1,i2,i3))->second;
  }
  virtual Outcome operator()(I1 i1, I2 i2, I3 i3) {
    return lookup(i1,i2,i3);
  }
};

/* RcppExport SEXP testTable(SEXP df, SEXP x) { */
/*   Table1D<double,double> table =  */
/*     Table1D<double,double>(as<DataFrame>(df)); */
/*   return wrap<double>(table.lookup(as<double>(x))); */
/* } */

/* RcppExport SEXP testTable2(SEXP df, SEXP x1, SEXP x2) { */
/*   Table2D<double,double,double> table =  */
/*     Table2D<double,double,double>(as<DataFrame>(df)); */
/*   return wrap<double>(table.lookup(as<double>(x1),as<double>(x2))); */
/* } */

/* RcppExport SEXP testTable3(SEXP df, SEXP x1, SEXP x2, SEXP x3) { */
/*   Table3D<double,double,double,double> table =  */
/*     Table3D<double,double,double,double>(as<DataFrame>(df)); */
/*   return wrap<double>(table.lookup(as<double>(x1),as<double>(x2), */
/* 				   as<double>(x3))); */
/* } */

#endif /* RCPP_TABLE_H */
