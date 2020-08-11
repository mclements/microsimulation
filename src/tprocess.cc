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
#include <siena/tprocess.h>

#if TPROCESS_IMPL==0
#else

#include <cassert>
#include <cstdio>

#if TPROCESS_IMPL==1
#include <ucontext.h>
#else
#include <setjmp.h>
#include <signal.h>
#endif

#include <siena/ssim.h>

namespace ssim {

static TProcess * current_tprocess = 0;

unsigned long TProcess::DefaultStackSize = 8 * 1024;

#if TPROCESS_IMPL==1
static ucontext_t resume_ctx;
#else
static jmp_buf resume_ctx; 
#endif

TProcess::TProcess(unsigned long size) {
    mystack_size = size;
    mystack = new char[mystack_size];
}

TProcess::TProcess() {
    mystack_size = DefaultStackSize;
    mystack = new char[mystack_size];
}

TProcess::~TProcess() {
    delete(mystack);
}

void TProcess::resume() {
#if TPROCESS_IMPL==1
    if (swapcontext(&resume_ctx, &running_ctx) < 0) {
	//
	// Error handling goes here
	//
	// ...work in progress...
	perror("TProcess::resume: swapcontext failed");
    }
#else
    if (!setjmp(resume_ctx)) {
	longjmp(running_ctx, 1);
    }
#endif
    //
    // returning from pause();
    //
}

void TProcess::pause() {
#if TPROCESS_IMPL==1
    if (swapcontext(&running_ctx, &resume_ctx) < 0) {
	//
	// Error handling goes here
	//
	// ...work in progress...
	perror("TProcess::pause: swapcontext failed");
	return;
    }
#else
    if (!setjmp(running_ctx)) {
	longjmp(resume_ctx, 1);
    }
#endif
    //
    // returning from resume();
    //
    current_tprocess = this;
}

#if TPROCESS_IMPL==1
void TProcess::starter() {
    TProcess * p = current_tprocess;
    p->main();
    Sim::stop_process();
    p->pause();
    // we should never get to this point, because the simulator should
    // never call process_event or process_timeout after we call
    // Sim::stop_process()
    assert(false);
}
#else
void TProcess::starter(int) {
    TProcess * p = current_tprocess;
    if (setjmp(p->running_ctx)) {
	p->main();
	Sim::stop_process();
	p->pause();
	// we should never get to this point, because the simulator should
	// never call process_event or process_timeout after we call
	// Sim::stop_process()
	assert(false);
    }
}
#endif

void TProcess::initialize() {
    //
    // this method creates the "execution context" for the thread that
    // executes this TProcess.  We can use two types of
    // implementations.  The first one is based on the "ucontext"
    // functions makecontext(), getcontext(), and swapcontext().  The
    // second implementation is based on a combination of POSIX
    // functions.
    // Breaking Change: this method was called init().
    //
#if TPROCESS_IMPL==1
    // this method uses makecontext and swapcontext and is pretty
    // straightforward.
    //
    if (getcontext(&running_ctx)) {	// first we create a running context
	perror("TProcess::initialize: getcontext failed for running_ctx");
	return;
    }
    running_ctx.uc_link = NULL;
    running_ctx.uc_stack.ss_sp = mystack;
    running_ctx.uc_stack.ss_size = mystack_size;
    running_ctx.uc_stack.ss_flags = 0;
    current_tprocess = this;
    makecontext(&running_ctx, TProcess::starter, 0);

    if (swapcontext(&resume_ctx, &running_ctx) < 0) {
	perror("TProcess::initialize: swapcontext failed");
	return;
    } 
#else
    // this method exploits a signal handler to create a different
    // execution context (on a different stack), and then uses
    // setjmp/longjmp to switch contexts.
    //
    stack_t stack_descr;
    struct sigaction sa;

    stack_descr.ss_flags = 0;		// set things up to use the new stack
    stack_descr.ss_size = mystack_size;
    stack_descr.ss_sp = mystack;
    sigaltstack(&stack_descr, 0);
    
    sa.sa_handler = TProcess::starter;	// set up the custom signal
    sa.sa_flags = SA_ONSTACK;		// handler to engage the starter
    sigemptyset(&sa.sa_mask); 
    sigaction(SIGUSR1, &sa, 0);
        
    current_tprocess = this;		// point the starter to this TProcess
    raise(SIGUSR1);			// object and fire up the starter

    resume();				// switch to the newly created context
#endif
    //
    // returning from pause();
    //
}

void TProcess::process_event(const Event * msg) {
    ev = msg;
    resume();
}

const Event * TProcess::wait_for_event(Time timeout) {
    Timeout * te = 0;
    if (timeout != INIT_TIME) {
	te = new Timeout();
	Sim::self_signal_event(te, timeout);
    }
    const Event * msg;
    do {
	current_tprocess->pause();
	msg = current_tprocess->ev;
	current_tprocess->ev = 0;
	//
	// we iterate to skip previously signaled Timeout
	//
    } while (dynamic_cast<const Timeout *>(msg) != 0 && msg != te);

    return msg;
}

void TProcess::stop() { };

}; // namespace ssim

#endif // TPROCESS_IMPL == 0
