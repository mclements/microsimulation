#ifndef RCPP_TABLE_H
#define RCPP_TABLE_H

#include <Rcpp.h>

#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <set>

using namespace std;
using namespace Rcpp;
using boost::get;

// TODO: re-write the DataFrameView class and adapt the Table class for it.

class Interpolate {
 public:
  vector<double> x, y, slope;
  Interpolate()  {
  }
 Interpolate(vector<double> inx, vector<double> iny) : 
  x(inx), y(iny) { 
    // calculate the slope between points
    for (size_t i=0; i<x.size()-1; i++) {
      slope.push_back((y[i+1]-y[i]) / (x[i+1]-x[i]));
    }
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

class NumericInterpolate {
 public:
  NumericVector x, y, slope;
  int n;
  NumericInterpolate()  {
  }
  NumericInterpolate(DataFrame df, int i0=0, int i1=1) { 
    // calculate the slope between points
    x = df(i0);
    y = df(i1);
    n = x.size();
    for (int i=0; i<n-1; i++) {
      slope.push_back((y[i+1]-y[i]) / (x[i+1]-x[i]));
    }
  }
  double approx(double xfind) {
    int i;
    if (xfind<=x[0]) return y[0];
    else if (xfind>=x[n-1]) return y[n-1];
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

template<class T>
T set_lower_bound(set<T,greater<T> > aset, T value) {
  return *aset.lower_bound(value);
}

template <class T>
class DataFrameSelect {
 public:
  Vector<Rcpp::traits::r_sexptype_traits<T>::rtype> data; 
  DataFrameSelect(const DataFrame & df, int i = 0) {
    data = df(i); // copy
  }
  DataFrameSelect(const DataFrame & df, string name) {
    data = df[name];
  }
  T operator[](int i) {
    return(data[i]);
  }
  int size() {
    return data.size();
  }
};

template<class Index0, class Outcome>
  class Table {
 public:
  typedef Index0 Index;
  typedef set<Index, greater<Index> > Axis;
  void insert(Index index, Outcome outcome) {
    axis.insert(index);
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[*axis.lower_bound(index)];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1) { 
    DataFrameSelect<Index0> df0(df,s0);
    DataFrameSelect<Outcome> df1(df,s1);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i]), df1[i]);
    }
  }
 private:
  Axis axis;
  map<Index,Outcome> data;
};

template<class Index0, class Index1, class Outcome>
  class Table<pair<Index0,Index1>,Outcome> {
 public:
  typedef boost::tuple<Index0,Index1> Index;
  typedef boost::tuple<
    set<Index0, greater<Index0> >,
    set<Index1, greater<Index1> > 
    > Axis;
  void insert(Index index, Outcome outcome) {
    get<0>(axis).insert(get<0>(index));
    get<1>(axis).insert(get<1>(index));
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[Index(*get<0>(axis).lower_bound(get<0>(index)),
		      *get<1>(axis).lower_bound(get<1>(index)))];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2) { 
    DataFrameSelect<Index0> df0(df,s0);
    DataFrameSelect<Index1> df1(df,s1);
    DataFrameSelect<Outcome> df2(df,s2);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i],df1[i]), df2[i]);
    }
  }
 private:
  Axis axis;
  map<Index,Outcome> data;
};

template<class I0, class I1, class Outcome>
  class Table<boost::tuple<I0,I1>,Outcome> {
 public:
  typedef boost::tuple<I0,I1> Index;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> > 
    > Axis;
  void insert(Index index, Outcome outcome) {
    get<0>(axis).insert(get<0>(index));
    get<1>(axis).insert(get<1>(index));
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[Index(*get<0>(axis).lower_bound(get<0>(index)),
		      *get<1>(axis).lower_bound(get<1>(index)))];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<Outcome> df2(df,s2);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i],df1[i]), df2[i]);
    }
  }
 private:
  Axis axis;
  map<Index,Outcome> data;
};


template<class I0, class I1, class I2, class Outcome>
  class Table<boost::tuple<I0,I1,I2>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2> Index;
  void insert(Index index, Outcome outcome) {
    axis0.insert(get<0>(index));
    axis1.insert(get<1>(index));
    axis2.insert(get<2>(index));
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[Index(*axis0.lower_bound(get<0>(index)),
		      *axis1.lower_bound(get<1>(index)),
		      *axis2.lower_bound(get<2>(index))
		      )];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2, string s3) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<Outcome> df3(df,s3);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i],df1[i],df2[i]), df3[i]);
    }
  }
 private:
  set<I0, greater<I0> > axis0;
  set<I1, greater<I1> > axis1;
  set<I2, greater<I2> > axis2;
  map<Index,Outcome> data;
};

template<class I0, class I1, class I2, class I3, class Outcome>
  class Table<boost::tuple<I0,I1,I2,I3>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2,I3> Index;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >
    > Axis;
  void insert(Index index, Outcome outcome) {
    get<0>(axis).insert(get<0>(index));
    get<1>(axis).insert(get<1>(index));
    get<2>(axis).insert(get<2>(index));
    get<3>(axis).insert(get<3>(index));
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[Index(*get<0>(axis).lower_bound(get<0>(index)),
		      *get<1>(axis).lower_bound(get<1>(index)),
		      *get<2>(axis).lower_bound(get<2>(index)),
		      *get<3>(axis).lower_bound(get<3>(index))
		      )];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2, string s3, string s4) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<I3> df3(df,s3);
    DataFrameSelect<Outcome> df4(df,s4);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i],df1[i],df2[i],df3[i]), df4[i]);
    }
  }
 private:
  Axis axis;
  map<Index,Outcome> data;
};


template<class I0, class I1, class I2, class I3, class I4, class Outcome>
  class Table<boost::tuple<I0,I1,I2,I3,I4>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2,I3,I4> Index;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >,
    set<I4, greater<I4> >
    > Axis;
  void insert(Index index, Outcome outcome) {
    get<0>(axis).insert(get<0>(index));
    get<1>(axis).insert(get<1>(index));
    get<2>(axis).insert(get<2>(index));
    get<3>(axis).insert(get<3>(index));
    get<4>(axis).insert(get<4>(index));
    data[index] = outcome;
  }
  virtual Outcome operator()(Index index) {
    return data[Index(*get<0>(axis).lower_bound(get<0>(index)),
		      *get<1>(axis).lower_bound(get<1>(index)),
		      *get<2>(axis).lower_bound(get<2>(index)),
		      *get<3>(axis).lower_bound(get<3>(index)),
		      *get<4>(axis).lower_bound(get<4>(index))
		      )];
  }
  Table(const DataFrame & df, string s0, string s1, string s2, string s3, string s4, string s5) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<I3> df3(df,s3);
    DataFrameSelect<I4> df4(df,s4);
    DataFrameSelect<Outcome> df5(df,s5);
    for (int i=0; i<df0.size(); i++) {
      insert(Index(df0[i],df1[i],df2[i],df3[i],df4[i]), df5[i]);
    }
  }
 private:
  Axis axis;
  map<Index,Outcome> data;
};

#endif /* RCPP_TABLE_H */
