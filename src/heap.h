// -*-C++-*-
//
//  This file is part of SSim, a simple discrete-event simulator.
//  See http://www.inf.unisi.ch/carzaniga/ssim/
//
//  Author: Antonio Carzaniga <firstname.lastname@unisi.ch>
//  See the file AUTHORS for full details. 
//
//  Copyright (C) 2004-2005 University of Colorado
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the GNU General Public License
//  as published by the Free Software Foundation; either version 2
//  of the License, or (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307,
//  USA, or send email to serl@cs.colorado.edu.
//
//
// $Id: heap.h,v 1.4 2005/10/16 21:01:49 carzanig Exp $
//
#ifndef _heap_h
#define _heap_h

#include <vector>

namespace ssim {

//
// This is an implementation of a binary heap taken from R. Sedgewick,
// "Algorithms in C," 3rd Ed., Vol. 1, pp. 368--375.  A few
// modifications to Sedgewick's implementation are noted below.
// Notice that this implementation is completely independent from the
// rest of SSim, and in fact can be used as is elsewhere.
//
template <typename T>
class heap {
public:
    typedef typename std::vector<T>::size_type		size_type;
    typedef typename std::vector<T>::iterator		iterator;
    typedef typename std::vector<T>::const_iterator	const_iterator;

private:
    std::vector<T> a;
    static void swap(T & first, T & second) throw() {
	T tmp = first;
	first = second;
	second = tmp;
    } 
    // contrary to Sedgewick's implementation, which counts elements
    // starting from 1, we start from 0.  So, in order to avoid
    // confusion and mistakes, we abstract all the position stuff with
    // the following member functions.  The compiler should optimize
    // everything out anyway.
    //
    static const size_type FIRST = 0;
    size_type last() const throw() { return a.size() - 1; }
    static size_type left(size_type pos) throw() { return pos*2 + 1; }
    static size_type right(size_type pos) throw() { return (pos + 1)*2; }
    static size_type parent(size_type pos) throw() { return (pos - 1)/2; }

public:
    bool empty() throw() { return a.empty(); }
    iterator begin() throw() { return a.begin(); }
    iterator end() throw() { return a.end(); }
    const_iterator begin() const  throw() { return a.begin(); }
    const_iterator end() const throw() { return a.end(); }
    void clear() throw() { a.clear(); }
    iterator erase(iterator first, iterator last) throw() {return a.erase(first, last); }

    void insert(const T & x) throw() {
	a.push_back(x);
	size_type k = last();
	size_type k_parent;
	while(k > FIRST) {	
	    k_parent = parent(k);
	    if (a[k] < a[k_parent]) {
		swap(a[k], a[k_parent]);
		k = k_parent;
	    } else {
		return;
	    }
	}
    }

    T pop_first() throw() {
	// ASSERT( FIRST <= last() ).  I.e., ! empty()
	T res = a[FIRST];
	if (FIRST == last()) {
	    a.pop_back();
	    return res;
	}
	a[FIRST] = a[last()];
	a.pop_back();
	size_type k = FIRST;
	size_type k_next;
	for(;;) {
	    k_next = left(k);
	    if (k_next > last()) {
		break;
	    }
	    if (right(k) <= last() && a[right(k)] < a[k_next]) {
		k_next = right(k);
	    }
	    if (a[k_next] < a[k]) {
		swap(a[k], a[k_next]);
		k = k_next;
	    } else {
		break;
	    }
	}
	return res;
    }
};

} // end namespace ssim

#endif /* _ssim_h */

