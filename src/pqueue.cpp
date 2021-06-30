#include <Rcpp.h>
#include<pqueue.h>

using Element = ssim::pqueueElement<SEXP>;
using Queue = ssim::pqueue<SEXP>;
using Rcpp::as;
using Rcpp::XPtr;
using Rcpp::wrap;

/**
   Create a new priority queue
*/
RcppExport
SEXP pqueue__new(SEXP _lower) {
  bool lower = as<bool>(_lower);
  XPtr<Queue> ptr(new Queue(lower));
  return wrap(ptr);
}
/**
   Push an active event onto the priority queue
*/
RcppExport
SEXP pqueue__push(SEXP _ptr, SEXP _priority, SEXP event) {
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  double priority = as<double>(_priority);
  ptr->push(priority, event);
  return R_NilValue;
}
/**
   Pop an active event from the priority queue.
*/
RcppExport
SEXP pqueue__pop(SEXP _ptr) {
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  Element element = std::move(ptr->pop());
  return Rcpp::List::create(Rcpp::_["priority"]=element.priority, Rcpp::_["event"]=element.event);
}
/**
   Cancel events that satisfy a predicate (that is, a test)
*/
RcppExport
SEXP pqueue__cancel(SEXP _ptr, SEXP _predicate) {
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  Rcpp::Function predicate = as<Rcpp::Function>(_predicate);
  if (!ptr->empty())
    for(std::vector<Element>::size_type i=0; i < ptr->_elements.size(); i++) {
      if (as<bool>(predicate(ptr->_elements[i].event)))
	ptr->_elements[i].active = false;
    }
  ptr->_anyCancelled = true;
  return R_NilValue;
}
/**
   Check whether a priority queue is empty
*/
RcppExport
SEXP pqueue__empty(SEXP _ptr) {
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  return wrap(ptr->empty());
}
/**
   Clear a priority queue
*/
RcppExport
SEXP pqueue__clear(SEXP _ptr) {
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  ptr->clear();
  return R_NilValue;
}
