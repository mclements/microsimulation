## Process this file with automake to produce Makefile.in

info_TEXINFOS = rngstreams.texi

rngstreams_TEXINFOS = \
 rngstreams.texi \
 version.texi

noinst_DATA = 

EXTRA_DIST = \
 rngstreams.html \
 rngstreams.txt \
 rngstreams.pdf \
 texinfo.tex 

# make plain text file ...
rngstreams.txt: rngstreams.texi version.texi
	$(MAKEINFO) $(AM_MAKEINFOFLAGS) --fill-column=80 --no-headers --number-sections \
	        rngstreams.texi -o rngstreams.txt

# abbreviations
.PHONY: txt 
txt: rngstreams.txt

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in \
 mdate-sh \
 rngstreams.html \
 rngstreams.txt
