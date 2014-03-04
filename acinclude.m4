dnl
dnl  This file is part of SSim, a simple discrete-event simulator.
dnl  See http://www.inf.usi.ch/carzaniga/ssim/
dnl  
dnl  Copyright (C) 1998-2005 University of Colorado
dnl  Copyright (C) 2012 Antonio Carzaniga
dnl  
dnl  Authors: Antonio Carzaniga <firstname.lastname@usi.ch>
dnl           See AUTHORS for full details.
dnl  
dnl  SSim is free software: you can redistribute it and/or modify it under
dnl  the terms of the GNU General Public License as published by the Free
dnl  Software Foundation, either version 3 of the License, or (at your
dnl  option) any later version.
dnl  
dnl  SSim is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl  GNU General Public License for more details.
dnl  
dnl  You should have received a copy of the GNU General Public License
dnl  along with SSim.  If not, see <http://www.gnu.org/licenses/>.
dnl 
dnl 
dnl AC_OPT_PROFILING
dnl
AC_DEFUN([AC_OPT_PROFILING], [
AC_ARG_ENABLE(profiling, 
  AC_HELP_STRING([--enable-profiling],
	[include profiling information. Values are "yes", "coverage" and "no" (default is "no")]),
dnl this is to optionally compile with profiling
dnl I don't know too much about this, but it looks like
dnl -pg only works with static libraries, so I'm going to 
dnl disable shared libraries here.
  [ case "$enableval" in
        coverage )
	    CFLAGS_prof='-pg -fprofile-arcs -ftest-coverage'
	    CXXFLAGS_prof='-pg -fprofile-arcs -ftest-coverage'
	    LDFLAGS_prof='-pg'
	    LIBS_prof='-lgcov'
	    AC_MSG_RESULT([Enabling profiling and coverage information])
	    ;;
        * )
	    CFLAGS_prof='-pg'
	    CXXFLAGS_prof='-pg'
	    LDFLAGS_prof='-pg'
	    LIBS_prof=''
	    AC_MSG_RESULT([Enabling profiling information])
	    ;;
    esac
    AC_DISABLE_SHARED ], 
  [ CFLAGS_prof=''
    CXXFLAGS_prof=''
    LDFLAGS_prof=''
    LIBS_prof=''
    AC_ENABLE_SHARED ])
AC_SUBST(CFLAGS_prof)
AC_SUBST(CXXFLAGS_prof)
AC_SUBST(LDFLAGS_prof)
AC_SUBST(LIBS_prof)
])
dnl
dnl AC_TPROCESS_IMPL
dnl
AC_DEFUN([AC_TPROCESS_IMPL], [
ac_tprocess_impl=ucontext
AC_ARG_ENABLE(tprocess, 
  AC_HELP_STRING([--enable-tprocess],
	[Enables the TProcess feature (default is to use ucontext)]),
  [ ac_tprocess_impl="$enableval"])
if test "$ac_tprocess_impl" = yes -o "$ac_tprocess_impl" = ucontext; then
  AC_CHECK_FUNCS(swapcontext, , [ac_tprocess_impl=longjmp])
fi
if test "$ac_tprocess_impl" = longjmp; then
  AC_CHECK_FUNCS(longjmp sigaction sigemptyset sigaltstack, , [ac_tprocess_impl=no])
fi
tprocess_impl=0
case "$ac_tprocess_impl" in
     no)
       tprocess_impl=0
       AC_DEFINE([TPROCESS_IMPL], [0], [Implements TProcess with ucontext])
       AC_MSG_WARN([disabling TProcess implementation])
       ;;
     yes | ucontext )
       tprocess_impl=1
       AC_DEFINE([TPROCESS_IMPL], [1], [Implements TProcess with ucontext])
       AC_MSG_NOTICE([using TProcess implementation based on ucontext.h])
       ;;
     longjmp)
       tprocess_impl=2
       AC_DEFINE([TPROCESS_IMPL], [2], [Implements TProcess with ucontext])
       AC_MSG_NOTICE([using TProcess implementation based on setjmp/longjmp])
       ;;
     *)
       AC_MSG_ERROR([unknown TProcess implementation: $ac_tprocess_impl (valid values: yes, ucontext, longjmp, no)])
       ;;
esac
AC_SUBST(tprocess_impl)
AM_CONDITIONAL(TPROCESS_TESTS, test "$ac_tprocess_impl" != no)
])

