## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libparser.la

libparser_la_SOURCES = \
 parser.c  parser.h  parser_source.h \
 stringparser.c stringparser_lists.ch \
 functparser.c \
 functparser_source.h \
 functparser_struct.h \
 functparser_symbols.h \
 functparser_stringgen.ch \
 functparser_debug.ch \
 functparser_deriv.ch \
 functparser_eval.ch \
 functparser_init.ch \
 functparser_parser.ch \
 functparser_scanner.ch

EXTRA_DIST = \
 make_stringparser.pl \
 functparser_doc.dh \
 stringparser_doc.dh \
 stringparser_lists.ch.in

parser_srcdir = $(top_srcdir)/src/parser

# make stringparser.c (process file stringparser.c.in)
sinclude .dep-stringparser_c
stringparser.c: stringparser_lists.ch
stringparser_lists.ch: make_stringparser.pl stringparser_lists.ch.in
	$(parser_srcdir)/make_stringparser.pl \
		< $(parser_srcdir)/stringparser_lists.ch.in \
		> $(parser_srcdir)/stringparser_lists.ch

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 stringparser_lists.ch \
 stringparser_doc.dh \
 .dep-stringparser_c \
 Makefile.in
