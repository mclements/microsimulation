#ifndef TRIPLE_H
#define TRIPLE_H

#include <utility>

using namespace std;

// http://stackoverflow.com/questions/8906545/how-to-initialize-a-vector-in-c
template <typename T, size_t N>
T* begin(T(&arr)[N]) { return &arr[0]; }
template <typename T, size_t N>
T* end(T(&arr)[N]) { return &arr[0]+N; }

template<class A, class B, class C>
class tuple3 : pair<A,pair<B,C> > { 
  // template <>
  // A get<0> () { return first; };
  // template <>
  // B get<1> () { return second.first; };
  // template <>
  // C get<2> () { return second.second; };
};

template<class A, class B, class C, class D>
class tuple4 : pair<A,pair<B,pair<C,D> > > { };

template<class A, class B, class C, class D, class E>
class tuple5 : pair<A,pair<B,pair<C,pair<D,E> > > > { };

/// triple holds three objects of arbitrary type.
template<class _T1, class _T2, class _T3>
    struct triple
    {
      typedef _T1 first_type;    ///<  @c first_type is the first bound type
      typedef _T2 second_type;   ///<  @c second_type is the second bound type
      typedef _T3 third_type;   ///<  @c third_type is the third bound type

      _T1 first;                 ///< @c first is a copy of the first object
      _T2 second;                ///< @c second is a copy of the second object
      _T3 third;                ///< @c third is a copy of the third object

      triple()
	: first(), second(), third() { }

      /** Three objects may be passed to a @c triple constructor to be copied.  */
      triple(const _T1& __a, const _T2& __b, const _T3& __c)
	: first(__a), second(__b), third(__c) { }

      /** There is also a templated copy ctor for the @c triple class itself.  */
      template<class _U1, class _U2, class _U3>
      triple(const triple<_U1, _U2, _U3>& __p)
        : first(__p.first),
          second(__p.second),
          third(__p.third) { }

    };

  /// Two triples of the same type are equal iff their members are equal.
template<class _T1, class _T2, class _T3>
    inline bool
operator==(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return __x.first == __y.first && __x.second == __y.second && 
	__x.third == __y.third; }

  /// <http://gcc.gnu.org/onlinedocs/libstdc++/manual/utilities.html>
template<class _T1, class _T2, class _T3>
    inline bool
    operator<(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return __x.first < __y.first
	|| (!(__y.first < __x.first) && __x.second < __y.second)
	|| (!(__y.first < __x.first) && !(__y.second < __x.second) && __x.third < __y.third); }

  /// Uses @c operator== to find the result.
  template<class _T1, class _T2, class _T3>
    inline bool
  operator!=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return !(__x == __y); }

  /// Uses @c operator< to find the result.
  template<class _T1, class _T2, class _T3>
    inline bool
    operator>(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return __y < __x; }

  /// Uses @c operator< to find the result.
  template<class _T1, class _T2, class _T3>
    inline bool
    operator<=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return !(__y < __x); }

  /// Uses @c operator< to find the result.
  template<class _T1, class _T2, class _T3>
    inline bool
    operator>=(const triple<_T1, _T2, _T3>& __x, const triple<_T1, _T2, _T3>& __y)
    { return !(__x < __y); }

  /**
   *  @brief A convenience wrapper for creating a triple from two objects.
   *  @param  x  The first object.
   *  @param  y  The second object.
   *  @return   A newly-constructed triple<> object of the appropriate type.
   *
   *  The standard requires that the objects be passed by reference-to-const,
   *  but LWG issue #181 says they should be passed by const value.  We follow
   *  the LWG by default.
   */
  // _GLIBCXX_RESOLVE_LIB_DEFECTS
  // 181.  make_triple() unintended behavior
  template<class _T1, class _T2, class _T3>
    inline triple<_T1, _T2, _T3>
  make_triple(_T1 __x, _T2 __y, _T3 __z)
  { return triple<_T1, _T2, _T3>(__x, __y, __z); }

#endif /* TRIPLE_H */
