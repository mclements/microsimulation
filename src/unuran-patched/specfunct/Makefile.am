## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libspecfunct.la

# Routines taken from
# Cephes Math Library Release 2.3:  June, 1995
# by Stephen L. Moshier:
# [ See also README.cephes ]
CEPHES = \
 cephes_source.h \
 gamma.c \
 igam.c \
 incbet.c \
 ndtr.c \
 ndtri.c \
 polevl.c

RMATH =

OTHER = \
 log1p.c

libspecfunct_la_SOURCES = $(CEPHES) $(RMATH) $(OTHER)

noinst_HEADERS = \
 unur_specfunct_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
