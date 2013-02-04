// -*-C++-*-
//
//  This file is part of SSim, a simple discrete-event simulator.
//  See http://www.inf.unisi.ch/carzaniga/ssim/
//
//  Author: Antonio Carzaniga <firstname.lastname@unisi.ch>
//  See the file AUTHORS for full details. 
//
//  Copyright (C) 1998-2004 University of Colorado
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
//  USA, or send email to the author.
//
//
// $Id: ssim.h,v 1.24 2005-12-13 09:50:03 carzanig Exp $
//
#ifndef _ssim_h
#define _ssim_h

// #include <boost/function.hpp>
#include "boost/boost/function.hpp"

/** \file ssim.h 
 *
 *  This header file defines the simulator API.
 **/

/** \namespace ssim
 *
 *  @brief name space for the Siena simulator.
 *
 *  This namespace groups all the types and functionalities associated
 *  with the Siena simulator.  These include:
 *
 *  <ol>
 *  <li>the simulator API
 *
 *  <li>the base classes for processes and events
 *
 *  <li>a few other utility classes
 *  </ol>
 **/
namespace ssim {

/** @brief version identifier for this ssim library
 **/
extern const char *	Version;

/** @brief process identifier type
 **/
typedef int		ProcessId;

/** @brief no process will be identified by NULL_PROCESSID
 **/
const ProcessId		NULL_PROCESSID = -1;

/** @brief virtual time type
 *
 *  This type represents the basic time in the virtual (simulated)
 *  world.  Being defined as an integer type, virtual time is a
 *  discrete quantity.  The actual semantics of the time unit is
 *  determined by the simulated application.  In other words, a time
 *  interval of 1 may be interpreted as one second, one year, or any
 *  other time interval, depending on the semantics of the simulated
 *  application.
 *
 *  @see Sim::advance_delay(Time),
 *       Sim::signal_event(ProcessId, const Event*, Time),
 *       and Sim::self_signal_event(const Event*, Time).
 **/
typedef double		Time;

/** @brief beginning of time
 **/
const Time		INIT_TIME = 0;


/** @brief basic event in the simulation.
 *
 *  This base class represents a piece of information or a signal
 *  exchanged between two processes through the simulator.  Every
 *  simulated event must inherit from this class.  For example:
 *
 *  \code
 *  class TextMessage : public Event {
 *  public:
 *      char * message;
 *
 *      TextMessage(const char * m) {
 *          message = strdup(m);
 *      }

 *      ~TextMessage() {
 *          free(message);
 *      }
 *  };
 *  \endcode
 *
 *  @see Sim::signal_event(),
 *       Sim::self_signal_event(),
 *       Process::process_event(const Event*),
 *       TProcess::wait_for_event(Time).
 **/
class Event {
 public:
			Event(): refcount(0) {};
    virtual		~Event() {};

 private:
    mutable unsigned refcount;
    friend class SimImpl;	// this is an opaque implementation class
    friend class Sim;		// these need to be friends to manage refcount
};

/** Type declaration for a predicate function
 *  used to test for particular events
 **/
  typedef boost::function<bool (const Event *)> EventPredicate;


/** @brief Virtual class (interface) representing processes running
 *  within the simulator.
 *
 *  A simulated process must implement this basic interface.
 **/
class Process {
 public:
    virtual		~Process() {};

    /** @brief action executed when the process is initialized.
     *
     *  This method is not a constructor.  It is rather an
     *  initialization step that is executed during the simulation
     *  when the process is created within the simulation through
     *  Sim::create_process.  This initialization is guaranteed to be
     *  executed before any event is signaled to the process.
     **/
    virtual void	init(void) {};

    /** @brief action executed in response to an event signaled to
     *  this process.
     *
     *  The Event object signaled through this method should not be
     *  used outside this method, other than by signaling it to other
     *  processes through Sim::signal_event().  In fact, the Event
     *  pointer passed to this method should be considered invalid
     *  outside this method, as the simulator proceeds to de-allocate
     *  every event object that is not signaled to any other process.
     *
     *  <p>The implementation of this method may specify the duration
     *  of the actions associated with this response using the \link
     *  Sim::advance_delay(Time) advance_delay\endlink method.  By
     *  default, the duration of an action is 0.
     *
     *  <p>The implementation of this method may use the C++
     *  dynamic_cast feature to select an action on the basis of the
     *  actual type of the Event object.  Here's an example that
     *  illustrates this feature:
     *
     *  \code
     *  class Payment : public Event {
     *  public:
     *      int dollars;
     *      //...
     *  };
     *
     *  class Delivery : public Event {
     *  public:
     *      int count;
     *      string model;
     *      //...
     *  };
     *
     *  class GuitarShop : public Process {
     *      //...
     *      virtual void process_event(const Event * e) {
     *          const Payment * p;
     *          const Delivery * d;
     *
     *          if ((p = dynamic_cast<const Payment *>(e)) != 0) { 
     *              cout << "I got a payment for $" << p->dollars << endl;
     *          } else if ((d = dynamic_cast<const Delivery *>(e)) != 0) {
     *              cout << "I got a delivery of " << d->count 
     *                   << " " << d->model << " guitars" << endl;
     *          } else {
     *              cerr << "GuitarShop: error: unknown event." << endl;
     *          }
     *      }
     *  };
     *  \endcode
     *
     *  @see Sim::advance_delay(Time).
     **/
    virtual void	process_event(const Event * msg) {};

    /** @brief executed when the process is explicitly stopped.
     *
     *  A process is stopped by a call to
     *  Sim::stop_process(ProcessId).  This method is executed
     *  immediately after the process has processed all the events
     *  scheduled before the call to Sim::stop_process.
     **/
    virtual void	stop(void) {};
};

/** @brief utility Process class providing a utility interface with the
 *  simulator.
 *
 *  This is a sligtly more advanced Process class that provides 
 *  automatic management of its own process id.
 **/
class ProcessWithPId : public Process {
 public:
    /** @brief activates this process within the simulator.
     *
     *  Creates a simulator process with this process object.  The
     *  same ProcessWithPId can be activated only once.
     *
     *  @return the ProcessId of the created simulator process. Or
     *  NULL_PROCESSID if this process object is not associated with a
     *  simulation process.
     *
     *  @see Sim::create_process(Process*)
     **/
    ProcessId		activate() throw();

    /** @brief process id of this process.
     *
     *  @return the id of the simulation process with this objectk,
     *  or NULL_PROCESSID if no process is associated with this
     *  object.
     **/
    ProcessId		pid() const throw();

    ProcessWithPId() throw();

 private:
    ProcessId process_id;
};

/** @brief an error handler for simulation errors.
 *
 *  Simulation errors occur when an event is scheduled for a process
 *  that is either terminated or busy processing other events.  These
 *  conditions may or may not represent application errors.  In any
 *  case, the simulator delegates the handling of these conditions to
 *  an error handler, which is implemented by extending this class.
 *
 *  @see Sim::set_error_handler(SimErrorHandler *)
 **/ 
class SimErrorHandler {
public:
    virtual ~SimErrorHandler() {}

    /** @brief handles a clear operation.
     *
     *  This method is called by Sim::clear().  This enables
     *  any counters or other internal state of the error 
     *  handler to be reset as necessary.
     *
     *  @see Sim::clear()
     **/
    virtual void clear() throw() {}

    /** @brief handles busy-process conditions.
     *
     *  A busy-process condition occurs when a process is scheduled to
     *  process an event at a time when it is still busy processing
     *  other events.
     *
     *  This method is executed within the simulation in the context
     *  of the busy process.  This means that the handler may access
     *  the simulation's \link Sim::clock() current time\endlink as
     *  well as all the interface functions of the simulator.  In
     *  particular, this method may also signal events to the busy
     *  process (e.g., with \link Sim::self_signal_event(const Event*)
     *  self_signal_event()\endlink).
     *
     *  @param p is the id of the busy process
     *  @param e is the scheduled event (possibly NULL)
     **/ 
    virtual void handle_busy(ProcessId p, const Event * e) throw() {}

    /** @brief handles terminated-process conditions.
     *
     *  A terminated-process condition occurs when a process is
     *  scheduled to process an event at a time when it has already
     *  terminated its execution.
     *
     *  This method is executed within the simulation in the context
     *  of the busy process.  This means that the handler may access
     *  the simulation's \link Sim::clock() current time\endlink as
     *  well as all the interface functions of the simulator.  In
     *  particular, this method may also signal events to the busy
     *  process (e.g., with \link Sim::self_signal_event(const Event*)
     *  self_signal_event()\endlink).
     *
     *  @param p is the id of the terminated process
     *  @param e is the scheduled event (possibly NULL)
     **/ 
    virtual void handle_terminated(ProcessId p, const Event * e) throw() {}
};

/** @brief a generic discrete-event sequential simulator
 *
 *  This class implements a generic discrete-event sequential
 *  simulator.  Sim maintains and executes a time-ordered schedule of
 *  actions (or discrete events).  
 *
 *  Notice that this class is designed to have only static members.
 *  It should therefore be seen and used more as a module than a
 *  class.  In practice, this means that Sim <em>does not define
 *  simulation objects</em>, but rather a single, static simulator.
 *  This design is based on practical considerations: a simulation
 *  object would allow one to handle multiple simulations within the
 *  same program at the same time, however that situation is neither
 *  common nor desirable.  Simulations tend to use a lot of memory and
 *  therefore having many of them is probably a bad idea.  Also,
 *  having multiple simulation objects would require every process
 *  object to maintain a reference to the corresponding simulation
 *  context (object), so that processes can access the virtual clock
 *  of their own simulation, their scheduler, etc.  This is both a
 *  waste of memory, for the many aditional references, and of CPU
 *  time for all the indirect calls to the simulation methods.
 **/ 
class Sim {
public:
    /** @brief creates a new process
     *
     *  Creates a new process with the given Process object.  This
     *  method schedules the execution of the \link Process::init()
     *  init\endlink method for the given object.
     *
     *  This method can be used safely within the simulation as well
     *  as outside the simulation.
     *
     *  @returns the process id of the new simulation process.
     *  @see Process::init()
     **/
    static ProcessId	create_process(Process *) throw();

    /** @brief stops the execution of a given process */
    static int		stop_process(ProcessId) throw();
    /** @brief stops the execution of the current process */
    static void		stop_process() throw();

   /** @brief clears out internal data structures
    *
    *  Resets the simulator making it available for a completely new
    *  simulation.  All scheduled actions are deleted together with
    *  the associated events.  All process identifiers returned by
    *  previoius invocations of \link create_process(Process*)
    *  create_process\endlink are invalidated by this operation.
    *  Notice however that it is the responsibility of the simulation
    *  programmer to delete process objects used in the simulation.
    **/
    static void		clear() throw();

    /** @brief signal an event to the current process immediately
     *
     *  Signal an event to \link this_process() this
     *  process\endlink.  The response is scheduled immediately (i.e.,
     *  at the \link Sim::clock() current time\endlink).
     *
     *  This method must be used within the simulation.  The effect of
     *  using this method outside the simulation is undefined.
     *
     *  @param e is the signaled event (possibly NULL)
     *
     *  @see signal_event(ProcessId, const Event *)
     *       and Process::process_event(const Event *).
     **/
    static void		self_signal_event(const Event * e) throw();

    /** @brief signal an event to the current process at the given time 
     *
     *  Signal a delayed event to the current process.  The response is
     *  scheduled with the given delay.
     *
     *  This method must be used within the simulation.  The effect of
     *  using this method outside the simulation is undefined.
     *
     *  @param e is the signaled event (possibly NULL)
     *
     *  @param delay is the delay from the \link Sim::clock() current
     *  time\endlink
     *
     *  @see signal_event() 
     *       and Process::process_event(const Event *).
     **/
    static void		self_signal_event(const Event * e, Time delay) throw();

    /** @brief signal an event to the given process immediately
     *
     *  Signal an event to the given process.  The response is
     *  scheduled immediately (i.e., at the \link Sim::clock() current
     *  time\endlink).
     *
     *  This method must be used within the simulation.  The effect of
     *  using this method outside the simulation is undefined.
     *
     *  @param p is the destination process. This must be a valid
     *  process id.  That is, a process id returned by 
     *  \link create_process(Process*) * create_process\endlink.  
     *  Using an invalid process id has undefined effects.
     *
     *  @param e is the signaled event (possibly NULL)
     *
     *  @see self_signal_event() 
     *       and Process::process_event(ProcessId, const Event *).
     **/
    static void		signal_event(ProcessId p, const Event * e) throw();

    /** @brief signal an event to the given process at the given time 
     *
     *  Signal a delayed event to the current process.  The response is
     *  scheduled with the given delay.
     *
     *  This method must be used within the simulation.  The effect of
     *  using this method outside the simulation is undefined.
     *
     *  @param p is the destination process. This must be a valid
     *  process id.  That is, a process id returned by 
     *  \link create_process(Process*) * create_process\endlink.  
     *  Using an invalid process id has undefined effects.
     *
     *  @param e is the signaled event (possibly NULL)
     *
     *  @param d is the signal delay starting from the \link
     *  Sim::clock() current time\endlink
     *
     *  @see self_signal_event() 
     *       and Process::process_event(const Event *).
     **/
    static void		signal_event(ProcessId p, const Event * e, Time d) throw();

    /** @brief advance the execution time of the current process.
     *
     *  This method can be used to specify the duration of certain
     *  actions, or certain steps within the same action.  For example:
     *
     *  \code
     *  class SomeProcess : public Process {
     *      //...
     *      virtual void process_event(const Event * e) {
     *          //
     *          // suppose this response is called at (virtual) time X
     *
     *          // ...do something here...
     *
     *          // the above actions have a default duration of 0,
     *          // therefore the following event is scheduled at time X + 5
     *          //
     *          signal_event(e, p, 5);
     *     
     *          advance_delay(10);
     *
     *          // ...do something else here...
     *
     *          // advance_delay(10) specifies a duration of 10, therefore 
     *          // the following (NULL) event is scheduled at time X + 15;
     *          //
     *          signal_event(NULL, 5);
     *      }
     *  };
     *  \endcode
     *
     *  Notice that a simulation process correspond to a single
     *  logical thread.  This means that Process::process_event() and
     *  TProcess::main() are intended to process one event at a time.
     *  Because of this chosen semantics, the use of
     *  advance_delay(Time) may result in the current process missing
     *  some events.  Referring to the above example, the overall
     *  duration of the execution step defined by
     *  SimpleProcess::process_event() is 15 time units.  Therefore, a
     *  SimpleProcess executing its process_event at time X would miss
     *  every event scheduled for it between time X and X + 15.  A
     *  semantically identical situation would occur for a sequential
     *  process.  For example:
     *
     *  \code
     *  class SomeTProcess : public TProcess {
     *      //...
     *      virtual void main() {
     *          // ...
     *          const Event * e;
     *          e = wait_for_event();
     *
     *          // now suppose wait_for_event() returns a signal at time X
     *          // and we do something with this signal...
     *          // and this "something" costs us 20 time units
     *          advance_delay(20);
     *
     *          // now, our virtual clock is X + 20, so we have missed 
     *          // all the signals between X and X + 20
     *          e = wait_for_signal();
     *      }
     *  };
     *  \endcode
     *
     *  In this example, the process misses all the events signaled
     *  within 20 time units after the first event.  This is because
     *  the process is busy working on the first event.
     *
     *  The application may program a handler for missed events by
     *  programming and registering a SimErrorHandler object.
     *
     *  This method must be used within the simulation.  The effect of
     *  using this method outside the simulation is undefined.
     *
     *  @see Sim::clock().
     *  @see SimErrorHandler
     **/
    static void		advance_delay(Time) throw();
    
    /** @brief returns the current process
     *
     *  Process id of the process that is currently scheduled by the
     *  simulation.  This method can be used by a Process in \link
     *  Process::process_event(const Event *e) process_event\endlink
     *  or by a TProcess in \link TProcess::main() main\endlink to
     *  figure out its own process id.
     *
     *  Example:
     *  \code
     *  class SimpleProcess : public Process {
     *      void SimpleProcess::process_event(const Event *) {
     *          cout << "My process id is: " << Sim::this_pocess() << endl;
     *      }
     *  };
     *  \endcode
     *
     *  @return process id of the current process, or \link
     *  ssim::NULL_PROCESSID NULL_PROCESSID\endlink if called outside
     *  the simulation.
     **/
    static ProcessId	this_process() throw();

    /** @brief returns the current virtual time for the current process
     *  
     *  Current virtual time.
     *  
     *  Example:
     *  \code
     *  void LoopProcess::process_event(const Event * e) {
     *      cout << "Here I am at time:" << Sim::clock() << endl;
     *      self_signal_event(NULL, 20)
     *      cout << "I just scheduled myself again for time: " 
     *           << Sim::clock() + 20 << endl;
     *  }
     *  \endcode
     *
     *  @return current virtual time for the current process.
     *  @see advance_delay(Time)
     **/
    static Time		clock() throw();
    
    /** @brief starts execution of the simulation */
    static void		run_simulation();
    /** @brief stops execution of the simulation */
    static void		stop_simulation() throw();

    /** @brief stops the execution of the simulation at the given time
     *
     *  This method sets the absolute (virtual) time at which the
     *  simulation will terminate.  This method can be used to limit
     *  the duration of a simulation even in the presence of
     *  schedulable actions.  When called with the (default)
     *  INIT_TIME, the simulation is set for normal termination, that
     *  is, the simulation terminates in the absence of schedulable
     *  actions.
     *
     *  @see stop_simulation()
     **/
    static void		set_stop_time(Time t = INIT_TIME) throw();

    /** @brief  registers a handler for simulation errors.
     *
     *  An error handler is a callback object that handles all
     *  simulation errors.
     *
     *  @see SimErrorHandler
     **/
    static void		set_error_handler(SimErrorHandler *) throw();

    static void remove_event(EventPredicate pred) throw();
};

} // end namespace ssim

#endif /* _ssim_h */

