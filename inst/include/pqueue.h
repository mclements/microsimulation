#ifndef PQUEUE_H
#define PQUEUE_H

#include <vector> // vector
#include <algorithm> // push_heap, pop_heap, make_heap, for_each
#include <exception> // exception
#include <stdexcept> // length_error
#include <utility> // move
#include <Rcpp.h> // Rf_error, forward_exception_to_r

namespace ssim {

  template<class Event>
  struct pqueueElement {
    double priority;
    long _order;
    bool active;
    Event event;
    pqueueElement(double priority, long order, bool active, Event event)
      : priority(priority), _order(order), active(active), event(event) { }
    virtual ~pqueueElement() = default;
  };

  template<class Event>
  struct pqueueElementComparator {
    bool smaller;
    pqueueElementComparator(bool smaller = true) : smaller(smaller) {}
    bool operator()(pqueueElement<Event> const& msg1, pqueueElement<Event> const& msg2) const {
      return (smaller ? msg1.priority > msg2.priority : msg1.priority < msg2.priority) ||
	(msg1.priority == msg2.priority && msg1._order > msg2._order);
    }
  };

  template<class Event>
  class pqueue {
  private:
    using Element = pqueueElement<Event>;
    using Comparator = pqueueElementComparator<Event>;
    using size_type = typename std::vector<Element>::size_type;
    Comparator _compare;
    long _entryOrder;
  public:
    std::vector<Element> _elements; // public for pqueue__cancel()
    bool _anyCancelled; // public for pqueue__cancel()
    explicit pqueue(bool smaller = true) {
      _compare = Comparator(smaller);
      _entryOrder = 0;
      _anyCancelled = false;
    }
    /**
       Push an elements to the priority queue
    */
    void push(Element element)
    {
      _elements.push_back(std::move(element));
      std::push_heap(_elements.begin(), _elements.end(), _compare);
    }
    void push(double priority, Event event)
    {
      Element element(priority, _entryOrder, true, event);
      push(std::move(element));
      // _elements.push_back(std::move(element));
      // std::push_heap(_elements.begin(), _elements.end(), _compare);
      _entryOrder++;
    }
    /**
       Pop an elements from the priority queue (that is, get the next element)
       Ignores whether an element is active or not.
    */
    Element pop()
    {
      try {
        if (empty()) {
          throw std::length_error("Empty priority queue");
        }
        while(!_elements.empty()) {
          std::pop_heap(_elements.begin(), _elements.end(), _compare);
          Element result = std::move(_elements.back());
          _elements.pop_back();
          if (result.active)
            return result; // change: std::move *not* used here
        }
      } catch(std::exception &ex) {
	forward_exception_to_r(ex);
      } catch(...) {
	::Rf_error("c++ exception (unknown reason)"); // R-specific
      }
      return _elements[0]; // never called (prevents -Wreturn-type warning)
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
       General method to apply a function f() to each element
       Use remake=true if f() may change the order in the priority queue.
    */
    template<class F>
    void for_each(F f, bool remake = false) {
      std::for_each(_elements.begin(), _elements.end(), f);
      if (remake)
	std::make_heap(_elements.begin(), _elements.end(), _compare);
    }
    /**
       Cancel elements that satisfy a predicate (that is, a test)
    */
    template<class Predicate>
    void cancel_element(Predicate predicate) {
      if (!empty()) {
	for(size_type i=0; i < _elements.size(); i++) {
	  if (predicate(_elements[i]))
	    _elements[i].active = false;
	}
      _anyCancelled = true;
      }
    }
    /**
       Cancel events that satisfy a predicate (that is, a test)
    */
    template<class Predicate>
    void cancel_event(Predicate predicate) {
      cancel_element([predicate](Element element) { return predicate(element.event); });
    }
    virtual ~pqueue() = default;
  };

  // SEXP pqueue__new(SEXP _lower);
  // SEXP pqueue__push(SEXP _ptr, SEXP _priority, SEXP event);
  // SEXP pqueue__pop(SEXP _ptr);
  // SEXP pqueue__cancel(SEXP _ptr, SEXP _predicate);
  // SEXP pqueue__empty(SEXP _ptr);
  // SEXP pqueue__clear(SEXP _ptr);

} // namespace ssim

#endif
