## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libuniform.la

libuniform_la_SOURCES = \
 urng_builtin.c urng_builtin.h mrg31k3p.c fish.c mstd.c \
 urng_fvoid.c urng_fvoid.h \
 urng_gsl.c urng_gsl.h \
 urng_gslqrng.c urng_gslqrng.h \
 urng_prng.c urng_prng.h \
 urng_randomshift.c urng_randomshift.h \
 urng_rngstreams.c urng_rngstreams.h

URNG_HEADER_FILES = \
 unuran_urng_gsl.h \
 unuran_urng_prng.h \
 unuran_urng_rngstreams.h

include_HEADERS = $(URNG_HEADER_FILES)

# make header files for installation
unuran_urng_gsl.h: urng_gsl.h urng_gslqrng.h
	$(top_srcdir)/scripts/make_urng_header.pl urng_gsl.h urng_gslqrng.h > unuran_urng_gsl.h

unuran_urng_prng.h: urng_prng.h
	$(top_srcdir)/scripts/make_urng_header.pl urng_prng.h > unuran_urng_prng.h

unuran_urng_rngstreams.h: urng_rngstreams.h
	$(top_srcdir)/scripts/make_urng_header.pl urng_rngstreams.h > unuran_urng_rngstreams.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 $(URNG_HEADER_FILES)
