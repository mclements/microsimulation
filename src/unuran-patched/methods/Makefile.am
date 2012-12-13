## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libmethods.la

# Files that contain deprecated routines
if ENABLE_DEPRECATED
DEPRECATED_FILES = \
 deprecated_methods.c deprecated_methods.h \
 deprecated_tdrgw.h \
 deprecated_vmt.c  deprecated_vmt.h  deprecated_vmt_struct.h
else
DEPRECATED_FILES = 
endif


libmethods_la_SOURCES = \
 $(DEPRECATED_FILES) \
\
 x_gen.c  x_gen.h  x_gen_struct.h x_gen_source.h \
\
 auto.c   auto.h   auto_struct.h \
 arou.c   arou.h   arou_struct.h \
 ars.c    ars.h    ars_struct.h \
 cext.c   cext.h   cext_struct.h \
 cstd.c   cstd.h   cstd_struct.h \
 dari.c   dari.h   dari_struct.h \
 dau.c    dau.h    dau_struct.h \
 dext.c   dext.h   dext_struct.h \
 dgt.c    dgt.h    dgt_struct.h \
 dsrou.c  dsrou.h  dsrou_struct.h \
 dss.c    dss.h    dss_struct.h \
 dstd.c   dstd.h   dstd_struct.h \
 empk.c   empk.h   empk_struct.h \
 empl.c   empl.h   empl_struct.h \
 gibbs.c  gibbs.h  gibbs_struct.h \
 hinv.c   hinv.h   hinv_struct.h \
 hist.c   hist.h   hist_struct.h \
 hitro.c  hitro.h  hitro_struct.h \
 hrb.c    hrb.h    hrb_struct.h \
 hrd.c    hrd.h    hrd_struct.h \
 hri.c    hri.h    hri_struct.h \
 itdr.c   itdr.h   itdr_struct.h \
 mcorr.c  mcorr.h  mcorr_struct.h \
 mixt.c   mixt.h   mixt_struct.h \
 mvstd.c  mvstd.h  mvstd_struct.h \
 mvtdr.c  mvtdr.h  mvtdr_struct.h \
          mvtdr_newset.ch mvtdr_init.ch mvtdr_sample.ch \
          mvtdr_info.ch mvtdr_debug.ch \
 ninv.c   ninv.h   ninv_struct.h \
          ninv_newset.ch ninv_init.ch ninv_sample.ch \
          ninv_newton.ch ninv_regula.ch \
          ninv_debug.ch ninv_info.ch \
 norta.c  norta.h  norta_struct.h \
 nrou.c   nrou.h   nrou_struct.h \
 pinv.c   pinv.h   pinv_struct.h \
          pinv_newset.ch pinv_init.ch pinv_sample.ch \
          pinv_prep.ch pinv_newton.ch \
          pinv_debug.ch pinv_info.ch \
 srou.c   srou.h   srou_struct.h \
 ssr.c    ssr.h    ssr_struct.h \
 tabl.c   tabl.h   tabl_struct.h tabl_newset.ch \
          tabl_init.ch tabl_sample.ch tabl_debug.ch tabl_info.ch \
 tdr.c    tdr.h    tdr_struct.h tdr_newset.ch \
          tdr_init.ch tdr_gw_init.ch tdr_ps_init.ch \
          tdr_sample.ch tdr_gw_sample.ch tdr_ps_sample.ch tdr_ia_sample.ch \
          tdr_debug.ch tdr_gw_debug.ch tdr_ps_debug.ch tdr_info.ch \
 unif.c   unif.h   unif_struct.h \
 utdr.c   utdr.h   utdr_struct.h \
 vempk.c  vempk.h  vempk_struct.h \
 vnrou.c  vnrou.h  vnrou_struct.h 

noinst_HEADERS = \
 unur_methods_source.h

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in

