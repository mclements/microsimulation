dnl Process this file with autoconf to produce a configure script.
dnl $Id: configure.ac 24 2006-09-07 07:06:56Z leydold $

dnl Initialize.
AC_PREREQ(2.59)
AC_INIT([RNGSTREAMS], [1.0.1], [unuran@statistik.wu-wien.ac.at])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/RngStream.h])
AC_CONFIG_HEADER([config.h])

dnl Library versioning (current:revision:age)
dnl See the libtool manual for an explanation of the numbers
dnl 
dnl Rngstreams-1.0  0:0:0

RNGSTREAMS_LT_VERSION="0:0:0"
AC_SUBST([RNGSTREAMS_LT_VERSION])

dnl add dmalloc support
AM_WITH_DMALLOC

dnl things required by automake
AC_PROG_MAKE_SET

dnl Check for which system.
AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_CC
AC_PROG_CPP
AC_PROG_INSTALL
AC_PROG_LN_S
AC_CHECK_TOOL(AR, ar, :)
# AC_CHECK_TOOL(RANLIB, ranlib, :)
AC_PROG_RANLIB
AC_DISABLE_SHARED
AC_PROG_LIBTOOL

dnl Checks for libraries.

dnl Checks for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([limits.h malloc.h sys/time.h])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE

dnl Checks for library functions.

dnl Set flags for compiler
if test X"$GCC" = Xyes ; then
        AC_SUBST(AM_CFLAGS,-Wall)
fi


AC_CONFIG_FILES([\
	Makefile \
	src/Makefile \
	doc/Makefile \
        tests/Makefile \
	examples/Makefile ])
AC_OUTPUT

dnl -----------------------------------------

