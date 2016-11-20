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

/**
    Class for numerical interpolation for x and y.
    Includes methods to read in x and y from a data-frame or from pairs of (x,y).
    Includes methods for linear approximation (approx, x->y) and inversion of increasing (invert)
    and decreasing (invert_decreasing) values (y->x).
    Includes an operator for a stepwise, left continuous function x->y.
 **/

class NumericInterpolate {
 public:
  NumericVector x, y, slope;
  int n;
 NumericInterpolate() : x(0), y(0), slope(0), n(0) {
  }
  NumericInterpolate(DataFrame df, int i0=0, int i1=1) {
    // calculate the slope between points
    x = df(i0);
    y = df(i1);
    n = x.size();
    prepare();
  }
  void prepare() {
    for (int i=0; i<n-1; i++) {
      slope.push_back((y[i+1]-y[i]) / (x[i+1]-x[i]));
    }
  }
  void push_back(pair<double,double> xy) {
    x.push_back(xy.first);
    y.push_back(xy.second);
    n++;
  }
  double approx(double xfind) {
    int i;
    if (xfind<=x[0]) return y[0];
    else if (xfind>=x[n-1]) return y[n-1]+slope[n-2]*(xfind-x[n-1]); // linear
    else {
      i = lower_bound(x.begin(), x.end(), xfind) - 1 - x.begin();
      return y[i]+slope[i]*(xfind-x[i]);
    }
  }
  double operator()(double xfind) {
    if (xfind<=x[0]) return y[0];
    int i = lower_bound(x.begin(), x.end(), xfind) - x.begin();
    return y[--i];
  }
  double invert(double yfind) { // assumes that the function is increasing
    int i;
    if (yfind<=y[0]) return x[0];
    else if (yfind>=y[n-1]) return x[n-1]+(yfind-y[n-1])/slope[n-2];
    else {
      i = lower_bound(y.begin(), y.end(), yfind) - 1 - y.begin();
      return x[i]+(yfind-y[i])/slope[i];
    }
  }
  double invert_decreasing(double yfind) { // assumes that the function is decreasing
    int i;
    if (yfind>=y[0]) return x[0];
    else if (yfind<y[n-1]) return x[n-1]+(yfind-y[n-1])/slope[n-2];
    else {
      i = lower_bound(y.begin(), y.end(), yfind, greater<double>()) - 1 - y.begin();
      return x[i]+(yfind-y[i])/slope[i];
    }
  }
};


template<class T>
T set_lower_bound(set<T,greater<T> > aset, T value) {
  return value<*aset.rbegin() ? *aset.rbegin() : *aset.lower_bound(value);
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

/** @brief A table class for lookups.  For the case of a single key,
    this is a small extension to std::map, including the ability to
    read columns from a DataFrame. Looking up a key which is less than the
    lowest key value will use the lowest key.

 **/
struct null_type {};

template<class I0 = null_type, class I1 = null_type, class I2 = null_type,
  class I3 = null_type, class I4 = null_type, class Outcome = null_type>
  class Table {
 public:
 Table() {}
  typedef boost::tuple<I0,I1,I2,I3,I4> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >,
    set<I4, greater<I4> >
    > Axis;
  void insert(I0 key0, I1 key1, I2 key2, I3 key3, I4 key4, Outcome outcome) {
    key_type key = key_type(key0,key1,key2,key3,key4);
    get<0>(axis).insert(key0);
    get<1>(axis).insert(key1);
    get<2>(axis).insert(key2);
    get<3>(axis).insert(key3);
    get<4>(axis).insert(key4);
    data[key] = outcome;
  }
  virtual Outcome operator()(I0 i0, I1 i1, I2 i2, I3 i3, I4 i4) {
    return data[key_type(set_lower_bound(get<0>(axis), i0),
			 set_lower_bound(get<1>(axis), i1),
			 set_lower_bound(get<2>(axis), i2),
			 set_lower_bound(get<3>(axis), i3),
			 set_lower_bound(get<4>(axis), i4))];
  }
  Table(const DataFrame & df, string s0, string s1, string s2, string s3, string s4, string s5) {
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<I3> df3(df,s3);
    DataFrameSelect<I4> df4(df,s4);
    DataFrameSelect<Outcome> df5(df,s5);
    for (int i=0; i<df0.size(); i++) {
      insert(df0[i],df1[i],df2[i],df3[i],df4[i], df5[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,mapped_type> data;
};

template<class I0, class I1, class I2, class I3, class Outcome>
  class Table<I0,I1,I2,I3,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2,I3> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >
    > Axis;
 Table() {}
  void insert(I0 key0, I1 key1, I2 key2, I3 key3, mapped_type outcome) {
    key_type key = key_type(key0,key1,key2,key3);
    get<0>(axis).insert(key0);
    get<1>(axis).insert(key1);
    get<2>(axis).insert(key2);
    get<3>(axis).insert(key3);
    data[key] = outcome;
  }
  virtual Outcome operator()(I0 i0, I1 i1, I2 i2, I3 i3) {
    return data[key_type(set_lower_bound(get<0>(axis), i0),
			 set_lower_bound(get<1>(axis), i1),
			 set_lower_bound(get<2>(axis), i2),
			 set_lower_bound(get<3>(axis), i3))];
  }
  Table(const DataFrame & df, string s0, string s1, string s2, string s3, string s4) {
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<I3> df3(df,s3);
    DataFrameSelect<Outcome> df4(df,s4);
    for (int i=0; i<df0.size(); i++) {
      insert(df0[i],df1[i],df2[i],df3[i], df4[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,mapped_type> data;
};

template<class I0, class I1, class I2, class Outcome>
  class Table<I0,I1,I2,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >
    > Axis;
 Table() {}
  void insert(I0 key0, I1 key1, I2 key2, mapped_type outcome) {
    key_type key = key_type(key0,key1,key2);
    get<0>(axis).insert(key0);
    get<1>(axis).insert(key1);
    get<2>(axis).insert(key2);
    data[key] = outcome;
  }
  virtual Outcome operator()(I0 i0, I1 i1, I2 i2) {
    return data[key_type(set_lower_bound(get<0>(axis), i0),
			 set_lower_bound(get<1>(axis), i1),
			 set_lower_bound(get<2>(axis), i2))];
  }
  Table(const DataFrame & df, string s0, string s1, string s2, string s3) {
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<Outcome> df3(df,s3);
    for (int i=0; i<df0.size(); i++) {
      insert(df0[i],df1[i],df2[i], df3[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,mapped_type> data;
};

template<class I0, class I1, class Outcome>
  class Table<I0,I1,Outcome> {
 public:
  typedef boost::tuple<I0,I1> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >
    > Axis;
 Table() {}
  void insert(I0 key0, I1 key1, mapped_type outcome) {
    key_type key = key_type(key0,key1);
    get<0>(axis).insert(key0);
    get<1>(axis).insert(key1);
    data[key] = outcome;
  }
  virtual Outcome operator()(I0 i0, I1 i1) {
    return data[key_type(set_lower_bound(get<0>(axis), i0),
			 set_lower_bound(get<1>(axis), i1))];
  }
  Table(const DataFrame & df, string s0, string s1, string s2) {
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<Outcome> df2(df,s2);
    for (int i=0; i<df0.size(); i++) {
      insert(df0[i],df1[i], df2[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,mapped_type> data;
};

template<class key_type, class mapped_type>
  class Table<key_type,mapped_type> {
 public:
  typedef set<key_type, greater<key_type> > Axis;
 Table() {}
  void insert(const key_type& key, const mapped_type& outcome) {
    axis.insert(key);
    data[key] = outcome;
  }
  virtual mapped_type operator()(key_type key) {
    return data[set_lower_bound(axis,key)];
  }
  Table(const DataFrame & df, string s0, string s1) {
    DataFrameSelect<key_type> df0(df,s0);
    DataFrameSelect<mapped_type> df1(df,s1);
    for (int i=0; i<df0.size(); i++) {
      insert(key_type(df0[i]), mapped_type(df1[i]));
    }
  }
 private:
  Axis axis;
  map<key_type,mapped_type> data;
};

#endif /* RCPP_TABLE_H */
