// -*-C++-*-
//
//  This file is part of Microsimulation package for R.
//  See http://github.com/mclements/microsimulation
//  
//  Authors: Mark Clements <firstname.lastname@ki.se>
//           See DESCRIPTION for full details.
//  
//  Microsimulation is free software: you can redistribute it and/or modify it under
//  the terms of the GNU General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your
//  option) any later version.
//  
//  Microsimulation is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with Microsimulation.  If not, see <http://www.gnu.org/licenses/>.
//
/** \mainpage Microsimulation C++ API Documentation

This documentation describes the C++ application programming interface
to <a href="http://github.com/mclements/microsimulation/">Microsimulation</a>, a
very simple discrete-event sequential simulation library for R and C++.  The C++ simulations build upon the to <a href="http://www.inf.usi.ch/carzaniga/ssim/">Ssim</a> simulation library.
The
simulator implemented by <em>ssim</em> executes <em>process
objects</em>.  Process objects can be programmed either as
<em>reactive</em> or <em>sequential</em> processes.  A reactive
process is programmed by a "callback" function that defines the
discrete execution steps of that process, performed in response to an
<em>event</em>. A sequential process is programmed as a traditional
sequential thread that can explicitly receive <em>events</em>.

<p>The events received by a process represent interactions with other
processes, activities scheduled by the process itself, or timeouts.
The simulation proceeds by scheduling the responses of each process to
the event signalled to that process.  During these execution steps, a
process may signal events to itself and to other processes,
immediately or with a delay, thereby scheduling other execution steps.
The simulation terminates when no more actions are scheduled.

<p>The ssim library consists of essentially two classes defined
within the \link ssim ssim\endlink namespace: \link ssim::Sim
Sim\endlink, which defines the interface to the simulator and \link
ssim::Process Process\endlink, which defines the interface and base
class for a reactive process.
User processes can be programmed by extending Process.

<p>\link ssim::Sim Sim\endlink offers the basic primitives for
signaling \link ssim::Event events\endlink, and for creating,
starting, and stopping processes.  \link ssim::Process Process\endlink
declares the execution steps scheduled when a process is \link
ssim::Process::initialize() started\endlink, \link
ssim::Process::process_event() signaled\endlink, and \link
ssim::Process::stop() stopped\endlink.  Notice that \link ssim::Sim
Sim\endlink defines a single, static simulation module, rather than a
class for simulation objects.  (See \link ssim::Sim Sim\endlink for
more detailed comments.)

<p>The execution of the simulation is based on a <em>virtual
clock</em> that represents the time in the simulated world.  The
virtual clock is simply a counter, therefore the time unit is
determined by the semantics of the simulated processes. The initial
value of the virtual clock is 0.  The passage of (virtual) time in the
simulated world is explicitly controlled by each process, essentially
in three ways:

<ul>

<li>by "sleeping".  That is, by scheduling a "timeout" event for
itself after a given interval.  The library does not provide an
explicit timeout event class, but rather it leaves that to the
application.  The easiest way to implement a timeout is to signal a
NULL event;

<li>by signaling events to other processes with a given delay (see
\link ssim::Sim::signal_event(ProcessId,const Event*,Time)
signal_event()\endlink);

<li>by explicitly declaring the duration of an action or an execution
step using the \link ssim::Sim::advance_delay(Time)
Sim::advance_delay(Time delay)\endlink method.

</ul>

<p>The documentation of \link ssim::Sim::advance_delay(Time)
advance_delay\endlink provides an in-depth discussion of the semantics
of the simulation in relation to virtual time.

<p>In addition to the basic Process class, the library provides a
utility class \link ssim::ProcessWithPId ProcessWithPId\endlink, that
automates some common procedures for process implementations.

*/

