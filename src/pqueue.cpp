#include <Rcpp.h>
#include<pqueue.h>

using namespace Rcpp;
using namespace ssim;

RcppExport
SEXP pqueue__new(SEXP _lower) {
  bool lower = as<bool>(_lower);
  XPtr<pqueue<SEXP> > ptr(new pqueue<SEXP>(lower));
  return wrap(ptr);
}
RcppExport
SEXP pqueue__push(SEXP _ptr, SEXP _priority, SEXP event) {
  XPtr<pqueue<SEXP> > ptr = as<XPtr<pqueue<SEXP> > >(_ptr);
  double priority = as<double>(_priority);
  ptr->push(priority, event);
  return R_NilValue;
}
/**
   Pop an active event from the priority queue.
*/
RcppExport
SEXP pqueue__pop(SEXP _ptr) {
  XPtr<pqueue<SEXP> > ptr = as<XPtr<pqueue<SEXP> > >(_ptr);
  pqueueElement<SEXP> element = std::move(ptr->pop());
  return List::create(_["priority"]=element.priority, _["event"]=element.event);
}
RcppExport
SEXP pqueue__cancel(SEXP _ptr, SEXP _predicate) {
  using Element = pqueueElement<SEXP>;
  using Queue = pqueue<SEXP>;
  XPtr<Queue> ptr = as<XPtr<Queue> >(_ptr);
  Function predicate = as<Function>(_predicate);
  if (!ptr->empty())
    for(std::vector<Element>::size_type i=0; i < ptr->_elements.size(); i++) {
      if (as<bool>(predicate(ptr->_elements[i].event)))
	ptr->_elements[i].active = false;
    }
  ptr->_anyCancelled = true;
  return R_NilValue;
}
RcppExport
SEXP pqueue__empty(SEXP _ptr) {
  XPtr<pqueue<SEXP> > ptr = as<XPtr<pqueue<SEXP> > >(_ptr);
  return wrap(ptr->empty());
}
RcppExport
SEXP pqueue__clear(SEXP _ptr) {
  XPtr<pqueue<SEXP> > ptr = as<XPtr<pqueue<SEXP> > >(_ptr);
  ptr->clear();
  return R_NilValue;
}
