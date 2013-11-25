#ifndef RCPP_TABLE_H
#define RCPP_TABLE_H

#include <Rcpp.h>

#include <functional>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

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

template <class T0, class T1>
  class DataFrameSelect<std::pair<T0,T1> > {
 public:
  typedef std::pair<T0,T1> Tuple;
  typedef std::pair<
    Vector<Rcpp::traits::r_sexptype_traits<T0>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T1>::rtype> > Data;
  Data data;
  DataFrameSelect(const DataFrame & df, int i0 = 0, int i1 = 1) {
    data = Data(df(i0),
		df(i1));
  }
  DataFrameSelect(const DataFrame & df, string name0, string name1) {
    data = Data(df(name0),
		df(name1));
  }
  Tuple operator[](int i) {
    return(Tuple((data.first)[i],
		 (data.second)[i]));
  }
};

template <class T0, class T1>
  class DataFrameSelect<boost::tuple<T0,T1> > {
 public:
  typedef boost::tuple<T0,T1> Tuple;
  typedef boost::tuple<
    Vector<Rcpp::traits::r_sexptype_traits<T0>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T1>::rtype> > Data;
  Data data;
  DataFrameSelect(const DataFrame & df, int i0 = 0, int i1 = 1) {
    data = Data(df(i0),
		df(i1));
  }
  DataFrameSelect(const DataFrame & df, string name0, string name1) {
    data = Data(df(name0),
		df(name1));
  }
  Tuple operator[](int i) {
    return(Tuple(get<0>(data)[i],
		 get<1>(data)[i]));
  }
};

template <class T0, class T1, class T2>
  class DataFrameSelect<boost::tuple<T0,T1,T2> > {
 public:
  typedef boost::tuple<T0,T1,T2> Tuple;
  typedef boost::tuple<
    Vector<Rcpp::traits::r_sexptype_traits<T0>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T1>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T2>::rtype> > Data;
  Data data;
  DataFrameSelect(const DataFrame & df, int i0 = 0, int i1 = 1, int i2 = 2) {
    data = Data(df(i0),
		df(i1),
		df(i2));
  }
  DataFrameSelect(const DataFrame & df, string name0, string name1, string name2) {
    data = Data(df(name0),
		df(name1),
		df(name2));
  }
  Tuple operator[](int i) {
    return(Tuple(get<0>(data)[i],
		 get<1>(data)[i],
		 get<2>(data)[i]));
  }
};

template <class T0, class T1, class T2, class T3>
  class DataFrameSelect<boost::tuple<T0,T1,T2,T3> > {
 public:
  typedef boost::tuple<T0,T1,T2,T3> Tuple;
  typedef boost::tuple<
    Vector<Rcpp::traits::r_sexptype_traits<T0>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T1>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T2>::rtype>,
    Vector<Rcpp::traits::r_sexptype_traits<T3>::rtype> > Data;
  Data data;
  DataFrameSelect(const DataFrame & df, int i0 = 0, int i1 = 1, int i2 = 2, int i3 = 3) {
    data = Data(df(i0),
		df(i1),
		df(i2),
		df(i3));
  }
  DataFrameSelect(const DataFrame & df, string name0, string name1, string name2,
		  string name3) {
    data = Data(df(name0),
		df(name1),
		df(name2),
		df(name3));
  }
  Tuple operator[](int i) {
    return(Tuple(get<0>(data)[i],
		 get<1>(data)[i],
		 get<2>(data)[i],
		 get<3>(data)[i]));
  }
};

template<class T>
int lookup_index(set<T, greater<T> > aset, T value) {
  typename set<T, greater<T> >::iterator it = aset.lower_bound(value);
  return it==aset.end() ? *--it : *it;
}

template<class Index, class Outcome>
  class Table {
 public:
  typedef map<Index,Outcome,greater<Index> > Map;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1) { 
    Vector<Rcpp::traits::r_sexptype_traits<Index>::rtype> df0(df,i0);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df1(df,i1);
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
  class Table<pair<I1,I2>,Outcome> {
public:
  typedef pair<I1,I2> Index;
  typedef map<Index,Outcome> Map;
  set<I1, greater<I1> > index1;
  set<I2, greater<I2> > index2;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1, int i2=2) { 
    DataFrameSelect<pair<I1,I2> > dfv(df,i0,i1);
    DataFrameSelect<Outcome> df2(df,i2);
    for (int i=0; i<df2.size(); i++) {
      Index index = dfv[i];
      index1.insert(index.first);
      index2.insert(index.second);
      data[index] = df2[i];
    }
  }
  virtual Outcome lookup(Index index) {
    Index new_index = Index(lookup_index(index1,index.first),
			    lookup_index(index2,index.second));
    return data.lower_bound(new_index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2) {
    return lookup(Index(i1,i2));
  }
  virtual Outcome operator()(I1 i1, I2 i2) {
    return lookup(Index(i1,i2));
  }
};

template<class I1, class I2, class Outcome>
  class Table<boost::tuple<I1,I2>,Outcome> {
public:
  typedef boost::tuple<I1,I2> Index;
  typedef map<Index,Outcome> Map;
  set<I1, greater<I1> > index1;
  set<I2, greater<I2> > index2;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1, int i2=2) { 
    DataFrameSelect<pair<I1,I2> > dfv(df,i0,i1);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df2(df,i2);
    for (int i=0; i<df2.size(); i++) {
      Index index = dfv[i];
      index1.insert(get<0>(index));
      index2.insert(get<1>(index));
      data[index] = df2[i];
    }
  }
  virtual Outcome lookup(Index index) {
    Index new_index = Index(lookup_index(index1,get<0>(index)),
			    lookup_index(index2,get<1>(index)));
    return data.lower_bound(new_index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2) {
    return lookup(Index(i1,i2));
  }
  virtual Outcome operator()(I1 i1, I2 i2) {
    return lookup(Index(i1,i2));
  }
};

template<class I1, class I2, class I3, class Outcome>
  class Table<boost::tuple<I1,I2,I3>,Outcome> {
public:
  typedef boost::tuple<I1,I2,I3> Index;
  typedef map<Index,Outcome,greater<Index> > Map;
  set<I1, greater<I1> > index1;
  set<I2, greater<I2> > index2;
  set<I3, greater<I3> > index3;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1, int i2=2, int i3=3) { 
    Vector<Rcpp::traits::r_sexptype_traits<I1>::rtype> df0 = df(i0);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df1 = df(i1);
    Vector<Rcpp::traits::r_sexptype_traits<I3>::rtype> df2 = df(i2);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df3 = df(i3);
    for (int i=0; i<df0.size(); i++) {
      Index index = Index(df0[i],df1[i],df2[i]);
      index1.insert(get<0>(index));
      index2.insert(get<1>(index));
      index3.insert(get<2>(index));
      data[index] = df3[i];
    }
  }
  virtual Outcome lookup(Index index) {
    Index new_index = Index(lookup_index(index1,get<0>(index)),
			    lookup_index(index2,get<1>(index)),
			    lookup_index(index3,get<2>(index)));
    return data.lower_bound(new_index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2, I3 i3) {
    return lookup(Index(i1,i2,i3));
  }
  virtual Outcome operator()(I1 i1, I2 i2, I3 i3) {
    return lookup(Index(i1,i2,i3));
  }
};

template<class I1, class I2, class I3, class I4, class Outcome>
  class Table<boost::tuple<I1,I2,I3,I4>,Outcome> {
public:
  typedef boost::tuple<I1,I2,I3,I4> Index;
  typedef map<Index,Outcome,greater<Index> > Map;
  set<I1, greater<I1> > index1;
  set<I2, greater<I2> > index2;
  set<I3, greater<I3> > index3;
  set<I4, greater<I4> > index4;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1, int i2=2, int i3=3, int i4=4) { 
    Vector<Rcpp::traits::r_sexptype_traits<I1>::rtype> df0 = df(i0);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df1 = df(i1);
    Vector<Rcpp::traits::r_sexptype_traits<I3>::rtype> df2 = df(i2);
    Vector<Rcpp::traits::r_sexptype_traits<I4>::rtype> df3 = df(i3);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df4 = df(i4);
    for (int i=0; i<df0.size(); i++) {
      Index index = Index(df0[i],df1[i],df2[i],df3[i]);
      index1.insert(get<0>(index));
      index2.insert(get<1>(index));
      index3.insert(get<2>(index));
      data[index] = df4[i];
    }
  }
  virtual Outcome lookup(Index index) {
    Index new_index = Index(lookup_index(index1,get<0>(index)),
			    lookup_index(index2,get<1>(index)),
			    lookup_index(index3,get<2>(index)),
			    lookup_index(index4,get<3>(index)));
    return data.lower_bound(new_index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2, I3 i3, I4 i4) {
    return lookup(Index(i1,i2,i3,i4));
  }
  virtual Outcome operator()(I1 i1, I2 i2, I3 i3, I4 i4) {
    return lookup(Index(i1,i2,i3,i4));
  }
};

template<class I1, class I2, class I3, class I4, class I5, class Outcome>
  class Table<boost::tuple<I1,I2,I3,I4,I5>,Outcome> {
public:
  typedef boost::tuple<I1,I2,I3,I4,I5> Index;
  typedef map<Index,Outcome,greater<Index> > Map;
  set<I1, greater<I1> > index1;
  set<I2, greater<I2> > index2;
  set<I3, greater<I3> > index3;
  set<I4, greater<I4> > index4;
  set<I5, greater<I5> > index5;
  Map data;
  Table(DataFrame df, int i0=0, int i1=1, int i2=2, int i3=3, int i4=4, int i5=5) { 
    Vector<Rcpp::traits::r_sexptype_traits<I1>::rtype> df0 = df(i0);
    Vector<Rcpp::traits::r_sexptype_traits<I2>::rtype> df1 = df(i1);
    Vector<Rcpp::traits::r_sexptype_traits<I3>::rtype> df2 = df(i2);
    Vector<Rcpp::traits::r_sexptype_traits<I4>::rtype> df3 = df(i3);
    Vector<Rcpp::traits::r_sexptype_traits<I5>::rtype> df4 = df(i4);
    Vector<Rcpp::traits::r_sexptype_traits<Outcome>::rtype> df5 = df(i5);
    for (int i=0; i<df0.size(); i++) {
      Index index = Index(df0[i],df1[i],df2[i],df3[i],df4[i]);
      index1.insert(get<0>(index));
      index2.insert(get<1>(index));
      index3.insert(get<2>(index));
      index4.insert(get<3>(index));
      index5.insert(get<4>(index));
      data[index] = df5[i];
    }
  }
  virtual Outcome lookup(Index index) {
    Index new_index = Index(lookup_index(index1,get<0>(index)),
			    lookup_index(index2,get<1>(index)),
			    lookup_index(index3,get<2>(index)),
			    lookup_index(index4,get<3>(index)),
			    lookup_index(index5,get<4>(index)));
    return data.lower_bound(new_index)->second;
  }
  virtual Outcome operator()(Index index) {
    return lookup(index);
  }
  virtual Outcome lookup(I1 i1, I2 i2, I3 i3, I4 i4, I5 i5) {
    return lookup(Index(i1,i2,i3,i4,i5));
  }
  virtual Outcome operator()(I1 i1, I2 i2, I3 i3, I4 i4, I5 i5) {
    return lookup(Index(i1,i2,i3,i4,i5));
  }
};

/* RcppExport SEXP testTable(SEXP df, SEXP x) { */
/*   Table1D<double,double> table =  */
/*     Table1D<double,double>(as<DataFrame>(df)); */
/*   return wrap<double>(table.lookup(as<double>(x))); */
/* } */

RcppExport SEXP testTable2(SEXP df, SEXP x1, SEXP x2) {
  Table<pair<double,double>,double> table =
    Table<pair<double,double>,double>(as<DataFrame>(df));
  return wrap<double>(table.lookup(as<double>(x1),as<double>(x2)));
}
// .Call("testTable2",data.frame(x1=rep(1:2,each=2),x2=rep(1:2,2),val=1:4),2.1,1.5,PACKAGE="microsimulation")


/* RcppExport SEXP testTable3(SEXP df, SEXP x1, SEXP x2, SEXP x3) { */
/*   Table3D<double,double,double,double> table =  */
/*     Table3D<double,double,double,double>(as<DataFrame>(df)); */
/*   return wrap<double>(table.lookup(as<double>(x1),as<double>(x2), */
/* 				   as<double>(x3))); */
/* } */

#endif /* RCPP_TABLE_H */
