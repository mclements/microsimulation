// -*-C++-*-
//
//  This file is part of SSim, a simple discrete-event simulator.
//  See http://www.inf.usi.ch/carzaniga/ssim/
//  
//  Copyright (C) 2003-2005 University of Colorado
//  Copyright (C) 2012 Antonio Carzaniga
//  
//  Authors: Antonio Carzaniga <firstname.lastname@usi.ch>
//           See AUTHORS for full details.
//  
//  SSim is free software: you can redistribute it and/or modify it under
//  the terms of the GNU General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your
//  option) any later version.
//  
//  SSim is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with SSim.  If not, see <http://www.gnu.org/licenses/>.
//
#ifndef _tprocess_h
#define _tprocess_h

#ifndef TPROCESS_IMPL
#include "tprocessconf.h"
#endif

// to make OSX happy with using ucontext.h
#define _XOPEN_SOURCE

#if TPROCESS_IMPL != 0

#if TPROCESS_IMPL==1
#include <ucontext.h>
#else
#include <setjmp.h>
#endif

#include <siena/ssim.h>

/** \file tprocess.h 
 *
 *  This header file defines \em sequential processes.
 **/

namespace ssim {

/** @brief Virtual class representing sequential processes.
 *
 *  A simulated \em sequential process extends this class by
 *  implementing a main() method.  The fundamental difference between
 *  the base Process and a TProcess is that a Process is implemented
 *  as a reactive object, that implements an algorithm by implementing
 *  reactions to incoming events in the form of callback functions.
 *  Conversely, a TProcess is a sequential process that implements an
 *  algorithm by simply implementing a main() function.
 *
 *  The \link tp.cc tp.cc example\endlink illustrates the use of
 *  sequential processes.
 *
 *  @see Process
 **/
class TProcess : public Process {
public:
    /** @brief creates a sequential process with the default stack size.
     * 
     *  @see DefaultStackSize
     **/
    TProcess();

    /** @brief creates a sequential process with the given stack size. */
    TProcess(unsigned long stacksize);

    virtual ~TProcess();

    /** @brief body of this simulated sequential process.
     *
     *  Body of this process.  This process may interact with other
     *  processes by \link wait_for_event(Time) receiving\endlink and
     *  \link Sim::signal_event() sending\endlink signals (a.k.a.,
     *  events).  The passage of time in this main body can be
     *  controlled through the \link Sim::advance_delay(Time)
     *  advance_delay\endlink method.
     **/
    virtual void main() = 0;

    /** @brief default stack size for sequential processes.
     *
     *  The initial value is 8Kb.  This value can be dynamically
     *  adjusted, and will affect all newly created TProcess objects.
     **/
    static unsigned long DefaultStackSize;

    /** @brief timeout event 
     *
     *  event type returned by wait_for_event(Time timeout) when the
     *  given timeout expires.
     **/
    class Timeout : public Event { };

    /** @brief accepts an event over a given interval.
     *
     *  This method suspends the execution of the current process
     *  until either an event is signaled to the current process or
     *  the given timeout expires.  A (default) timeout value of \link
     *  ssim::INIT_TIME INIT_TIME\endlink means an infinite timeout.
     *  
     *  <p>If an event is signaled before the timeout expires, this
     *  method returns that event.  If the timeout expires before any
     *  event is received, this method returns a Timeout event.
     *
     *  This method \em must be called within the execution of the
     *  main() method of a TProcess.  The following example
     *  illustrates a typical way to use this method:
     *
     *  \code
     *  class HaveFork : public Event {
     *      //...
     *  };
     *  
     *  class Philosopher : public TProcess {
     *      //...
     *      virtual void main() {
     *          cout << "I'll see if a fork becomes available" << endl;
     *          //
     *          const Event * e = wait_for_event(100);
     *          if (dynamic_cast<const HaveFork *>(e) != 0) {
     *              cout << "I got one, I'll start eating now..." << endl;
     *              // ...
     *          } else if (dynamic_cast<const TProcess::Timeout *>(e) != 0) {
     *              cout << "I got tired of waiting!" << endl;
     *              // ...
     *              cout << "I'll go to sleep now" << endl;
     *              cout << "Somebody must send me a signal to wake me up" << endl; 
     *              wait_for_event();
     *          } else {
     *              cout << "I'm getting confusing signals here..." << endl;
     *          }
     *      }
     *  };
     *  \endcode
     **/
    static const Event * wait_for_event(Time timeout = INIT_TIME);

private:
    virtual void init(void);
    virtual void process_event(const Event * msg);
    virtual void stop(void);

#if TPROCESS_IMPL==1
    ucontext_t running_ctx;
#else
    jmp_buf running_ctx; 
#endif
    char * mystack;
    unsigned long mystack_size;

    const Event * ev;

#if TPROCESS_IMPL==1
    static void starter();
#else
    static void starter(int);
#endif
    void pause();
    void resume();
};

}; // end namespace ssim

#endif /* TPROCESS_IMPL!=VOID */

#endif /* _ssim_h */

