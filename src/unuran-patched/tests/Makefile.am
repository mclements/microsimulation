## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libtests.la

libtests_la_SOURCES = \
 chi2test.c \
 correlation.c\
 countpdf.c \
 counturn.c \
 inverror.c \
 moments.c \
 printsample.c \
 quantiles.c \
 tests.c \
 timing.c

include_HEADERS = \
 unuran_tests.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
