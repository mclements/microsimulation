## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libdistr.la

# Files that contain deprecated routines
if ENABLE_DEPRECATED
DEPRECATED_FILES = deprecated_distr.c deprecated_distr.h 
else
DEPRECATED_FILES = 
endif

libdistr_la_SOURCES = \
 $(DEPRECATED_FILES) \
\
 distr.c   distr.h  distr_info.c distr_struct.h distr_source.h \
 cont.c    cont.h \
 corder.c  corder.h \
 cxtrans.c cxtrans.h \
 cemp.c    cemp.h \
 cvec.c    cvec.h \
 condi.c   condi.h \
 cvemp.c   cvemp.h\
 discr.c   discr.h \
 matr.c	   matr.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

