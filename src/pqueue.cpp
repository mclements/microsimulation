#include<Rcpp.h>
namespace ssim {
  
  using namespace Rcpp;
  
  struct pqueueElement {
    double priority;
    long _order;
    bool active;
    SEXP event;
    pqueueElement(double priority, long order, bool active, SEXP event)
      : priority(priority), _order(order), active(active), event(event) { }
    virtual ~pqueueElement() = default;
  };
  
  struct pqueueElementComparator {
    bool smaller;
    pqueueElementComparator(bool smaller = true) : smaller(smaller) {}
    bool operator()(pqueueElement const& msg1, pqueueElement const& msg2) const {
      return (smaller ? msg1.priority > msg2.priority : msg1.priority < msg2.priority) ||
	(msg1.priority == msg2.priority && msg1._order > msg2._order);
    }
  };
  
  class pqueue {
  private:
    std::vector<pqueueElement> _elements;
    pqueueElementComparator _compare;
    long _entryOrder;
    bool _anyCancelled;
    using size_type = std::vector<pqueueElement>::size_type;
  public:
    explicit pqueue(bool smaller = true) {
      _compare = pqueueElementComparator(smaller);
      _entryOrder = 0;
      _anyCancelled = false;
    }
    /**
       Push an elements to the priority queue
    */
    void pushElement(pqueueElement element)
    {
      _elements.push_back(std::move(element));
      std::push_heap(_elements.begin(), _elements.end(), _compare);
    }
    void push(double priority, SEXP event)
    {
      pqueueElement element(priority, _entryOrder, true, event);
      pushElement(std::move(element));
      // _elements.push_back(std::move(element));
      // std::push_heap(_elements.begin(), _elements.end(), _compare);
      _entryOrder++;
    }
    /**
       Pop an elements from the priority queue (that is, get the next element)
       Ignores whether an element is active or not.
    */
    pqueueElement popElement()
    {
      try {
        if (empty()) {
          throw std::length_error("Empty priority queue");
        }
        while(!_elements.empty()) {
          std::pop_heap(_elements.begin(), _elements.end(), _compare);
          pqueueElement result = std::move(_elements.back());
          _elements.pop_back();
          if (result.active) 
            return result; // change: std::move *not* used here
        }
      } catch(std::exception &ex) {	
	forward_exception_to_r(ex);
      } catch(...) { 
	::Rf_error("c++ exception (unknown reason)"); 
      }
      return _elements[0]; // never called (prevents -Wreturn-type warning)
    }
    /**
       Pop an active event from the priority queue.
    */
    List pop() {
      pqueueElement element = std::move(popElement());
      return List::create(_["priority"]=element.priority, _["event"]=element.event);
    }
    /**
       Check whether the priority queue is either empty or has only inactive events
    */
    bool empty() {
      if (_elements.empty()) return true; // empty queue
      else if (!_anyCancelled) return false; // queue not empty and no cancelled events
      else { // loop through to see if any of the elements are active
        for(size_type i=0; i<_elements.size(); i++)
          if(_elements[i].active) return false; // at least one active events
        return true; // no active events
      }
    }
    /**
       Clear the priority queue
    */
    void clear () {
      _elements.clear(); // assumes that elements are safe pointers
    }
    /**
       Cancel any events that satisfy a predicate (which is an R function).
    */
    void cancel(Rcpp::Function predicate) {
      if (!empty())
        for(size_type i=0; i<_elements.size(); i++) {
          if (as<bool>(predicate(_elements[i].event)))
            _elements[i].active = false;
        }
      _anyCancelled = true;
    }
    /**
       General method to apply a function f() to each element
       Use remake=true if f() may change the order in the priority queue.
    */
    template<class F>
    void for_each(F f, bool remake = false) {
      std::for_each(_elements.begin(), _elements.end(), f);
      if (remake)
	std::make_heap(_elements.begin(), _elements.end(), _compare);
    }
    virtual ~pqueue() = default;
  };
  
  RcppExport
  SEXP pqueue__new(SEXP _lower) {
    bool lower = as<bool>(_lower);
    Rcpp::XPtr<pqueue> ptr(new pqueue(lower));
    return wrap(ptr);
  }
  RcppExport
  SEXP pqueue__push(SEXP _ptr, SEXP _priority, SEXP event) {
    XPtr<pqueue> ptr = as<XPtr<pqueue> >(_ptr);
    double priority = as<double>(_priority);
    ptr->push(priority, event);
    return R_NilValue;
  }
  RcppExport
  SEXP pqueue__pop(SEXP _ptr) {
    XPtr<pqueue> ptr = as<XPtr<pqueue> >(_ptr);
    return wrap(ptr->pop());
  }
  RcppExport
  SEXP pqueue__cancel(SEXP _ptr, SEXP _predicate) {
    XPtr<pqueue> ptr = as<XPtr<pqueue> >(_ptr);
    Function predicate = as<Function>(_predicate);
    ptr->cancel(predicate);
    return R_NilValue;
  }
  RcppExport
  SEXP pqueue__empty(SEXP _ptr) {
    XPtr<pqueue> ptr = as<XPtr<pqueue> >(_ptr);
    return wrap(ptr->empty());
  }
  RcppExport
  SEXP pqueue__clear(SEXP _ptr) {
    XPtr<pqueue> ptr = as<XPtr<pqueue> >(_ptr);
    ptr->clear();
    return R_NilValue;
  }

} // namespace ssim
