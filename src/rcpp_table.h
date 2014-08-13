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
  double invert(double yfind) {
    int i;
    if (yfind<=y[0]) return x[0];
    else if (yfind>=y[n-1]) return x[n-1]+(yfind-y[n-1])/slope[n-2];
    else {
      i = lower_bound(y.begin(), y.end(), yfind) - 1 - y.begin();
      return x[i]+(yfind-y[i])/slope[i];
    }
  }
  double invert_decreasing(double yfind) {
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

template<class Key0, class Outcome>
  class Table {
 public:
  typedef Key0 key_type;
  typedef pair<key_type,Outcome> value_type;
  typedef set<key_type, greater<key_type> > Axis;
  typedef Outcome mapped_type;
  void insert(const key_type& key, const Outcome& outcome) {
    axis.insert(key);
    data[key] = outcome;
  }
  void push_back(const value_type& value) {
    axis.insert(value.first);
    data[value.first].push_back(value.second);
  }
  virtual Outcome operator()(key_type key) {
    return data[*axis.lower_bound(key)];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1) { 
    DataFrameSelect<Key0> df0(df,s0);
    DataFrameSelect<Outcome> df1(df,s1);
    for (int i=0; i<df0.size(); i++) {
      insert(key_type(df0[i]), df1[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,Outcome> data;
};

template<class Key0, class Key1, class Outcome>
  class Table<pair<Key0,Key1>,Outcome> {
 public:
  typedef boost::tuple<Key0,Key1> key_type;
  typedef Outcome mapped_type;
  typedef map<key_type,mapped_type> data_type;
  typedef boost::tuple<
    set<Key0, greater<Key0> >,
    set<Key1, greater<Key1> > 
    > Axis;
  void insert(key_type key, Outcome outcome) {
    get<0>(axis).insert(get<0>(key));
    get<1>(axis).insert(get<1>(key));
    data[key] = outcome;
  }
  virtual Outcome operator()(key_type key) {
    return data[key_type(*get<0>(axis).lower_bound(get<0>(key)),
		      *get<1>(axis).lower_bound(get<1>(key)))];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2) { 
    DataFrameSelect<Key0> df0(df,s0);
    DataFrameSelect<Key1> df1(df,s1);
    DataFrameSelect<Outcome> df2(df,s2);
    for (int i=0; i<df0.size(); i++) {
      insert(key_type(df0[i],df1[i]), df2[i]);
    }
  }
 private:
  map<key_type,mapped_type> data;
  Axis axis;
};

template<class I0, class I1, class Outcome>
  class Table<boost::tuple<I0,I1>,Outcome> {
 public:
  typedef boost::tuple<I0,I1> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> > 
    > Axis;
  void insert(key_type key, Outcome outcome) {
    get<0>(axis).insert(get<0>(key));
    get<1>(axis).insert(get<1>(key));
    data[key] = outcome;
  }
  virtual Outcome operator()(key_type key) {
    return data[key_type(*get<0>(axis).lower_bound(get<0>(key)),
			 *get<1>(axis).lower_bound(get<1>(key)))];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<Outcome> df2(df,s2);
    for (int i=0; i<df0.size(); i++) {
      insert(key_type(df0[i],df1[i]), df2[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,Outcome> data;
};


template<class I0, class I1, class I2, class Outcome>
  class Table<boost::tuple<I0,I1,I2>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2> key_type;
  typedef Outcome mapped_type;
  void insert(key_type key, Outcome outcome) {
    axis0.insert(get<0>(key));
    axis1.insert(get<1>(key));
    axis2.insert(get<2>(key));
    data[key] = outcome;
  }
  virtual Outcome operator()(key_type key) {
    return data[key_type(*axis0.lower_bound(get<0>(key)),
			 *axis1.lower_bound(get<1>(key)),
			 *axis2.lower_bound(get<2>(key))
			 )];
  }
  Table() {}
  Table(const DataFrame & df, string s0, string s1, string s2, string s3) { 
    DataFrameSelect<I0> df0(df,s0);
    DataFrameSelect<I1> df1(df,s1);
    DataFrameSelect<I2> df2(df,s2);
    DataFrameSelect<Outcome> df3(df,s3);
    for (int i=0; i<df0.size(); i++) {
      insert(key_type(df0[i],df1[i],df2[i]), df3[i]);
    }
  }
  map<key_type,Outcome> data;
 private:
  set<I0, greater<I0> > axis0;
  set<I1, greater<I1> > axis1;
  set<I2, greater<I2> > axis2;
};

template<class I0, class I1, class I2, class I3, class Outcome>
  class Table<boost::tuple<I0,I1,I2,I3>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2,I3> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >
    > Axis;
  void insert(key_type key, Outcome outcome) {
    get<0>(axis).insert(get<0>(key));
    get<1>(axis).insert(get<1>(key));
    get<2>(axis).insert(get<2>(key));
    get<3>(axis).insert(get<3>(key));
    data[key] = outcome;
  }
  virtual Outcome operator()(key_type key) {
    return data[key_type(*get<0>(axis).lower_bound(get<0>(key)),
		      *get<1>(axis).lower_bound(get<1>(key)),
		      *get<2>(axis).lower_bound(get<2>(key)),
		      *get<3>(axis).lower_bound(get<3>(key))
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
      insert(key_type(df0[i],df1[i],df2[i],df3[i]), df4[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,Outcome> data;
};


template<class I0, class I1, class I2, class I3, class I4, class Outcome>
  class Table<boost::tuple<I0,I1,I2,I3,I4>,Outcome> {
 public:
  typedef boost::tuple<I0,I1,I2,I3,I4> key_type;
  typedef Outcome mapped_type;
  typedef boost::tuple<
    set<I0, greater<I0> >,
    set<I1, greater<I1> >,
    set<I2, greater<I2> >,
    set<I3, greater<I3> >,
    set<I4, greater<I4> >
    > Axis;
  void insert(key_type key, Outcome outcome) {
    get<0>(axis).insert(get<0>(key));
    get<1>(axis).insert(get<1>(key));
    get<2>(axis).insert(get<2>(key));
    get<3>(axis).insert(get<3>(key));
    get<4>(axis).insert(get<4>(key));
    data[key] = outcome;
  }
  virtual Outcome operator()(key_type key) {
    return data[key_type(*get<0>(axis).lower_bound(get<0>(key)),
		      *get<1>(axis).lower_bound(get<1>(key)),
		      *get<2>(axis).lower_bound(get<2>(key)),
		      *get<3>(axis).lower_bound(get<3>(key)),
		      *get<4>(axis).lower_bound(get<4>(key))
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
      insert(key_type(df0[i],df1[i],df2[i],df3[i],df4[i]), df5[i]);
    }
  }
 private:
  Axis axis;
  map<key_type,Outcome> data;
};

#endif /* RCPP_TABLE_H */
