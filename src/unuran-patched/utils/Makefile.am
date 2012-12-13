## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libutils.la

libutils_la_DEPENDENCIES = unur_fp_const_source.h

libutils_la_SOURCES = \
 debug.c  debug.h debug_source.h \
 eigensystem.c \
 error.c  error.h error_source.h \
 fmax.c 	  fmax_source.h \
 hooke.c 	  hooke_source.h \
 lobatto.c        lobatto_source.h lobatto_struct.h \
 matrix.c 	  matrix_source.h \
 mrou_rectangle.c mrou_rectangle_source.h mrou_rectangle_struct.h \
 slist.c  slist.h slist_struct.h \
 stream.c stream.h stream_source.h \
 string.c         string_source.h string_struct.h \
 umalloc.c        umalloc_source.h \
 umath.c  umath.h umath_source.h \
 unur_fp.c  \
 vector.c 	  vector_source.h 

noinst_HEADERS = \
 unur_fp_source.h unur_fp_const_source.h \
 unur_errno.h \
 unur_math_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
##DISTCLEANFILES = \
## unur_fp_const_source.h

MAINTAINERCLEANFILES = \
 Makefile.in \
 unur_fp_const_source.h

