## Process this file with automake to produce Makefile.in

AM_LDFLAGS = -L$(top_builddir)/src

INCLUDES = -I$(top_srcdir)/src

LDADD = $(top_builddir)/src/librngstreams.la

EXAMPLES = \
 example1

noinst_PROGRAMS = $(EXAMPLES)

TESTS = $(EXAMPLES)

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = Makefile.in

