## Process this file with automake to produce Makefile.in
## $Id: Makefile.am 5330 2011-04-19 10:50:18Z leydold $

AM_CPPFLAGS = \
 -I$(top_srcdir)/src

noinst_LTLIBRARIES = libdistributions.la

libdistributions_la_SOURCES = \
 c_beta.c        c_beta_gen.c \
 c_burr.c \
 c_cauchy.c \
 c_chi.c         c_chi_gen.c \
 c_chisquare.c \
 c_exponential.c c_exponential_gen.c \
 c_extremeI.c \
 c_extremeII.c \
 c_F.c \
 c_gamma.c       c_gamma_gen.c \
 c_ghyp.c \
 c_gig.c         c_gig_gen.c \
 c_gig2.c \
 c_hyperbolic.c \
 c_ig.c \
 c_laplace.c \
 c_logistic.c \
 c_lognormal.c \
 c_lomax.c \
 c_normal.c      c_normal_gen.c \
 c_pareto.c \
 c_powerexponential.c c_powerexponential_gen.c \
 c_rayleigh.c \
 c_slash.c       c_slash_gen.c \
 c_student.c     c_student_gen.c \
 c_triangular.c \
 c_uniform.c \
 c_weibull.c \
 d_binomial.c    d_binomial_gen.c \
 d_geometric.c \
 d_hypergeometric.c d_hypergeometric_gen.c \
 d_logarithmic.c d_logarithmic_gen.c \
 d_negativebinomial.c \
 d_poisson.c     d_poisson_gen.c \
 d_zipf.c        d_zipf_gen.c \
 m_correlation.c \
 vc_copula.c \
 vc_multicauchy.c \
 vc_multiexponential.c \
 vc_multinormal.c   vc_multinormal_gen.c \
 vc_multistudent.c

noinst_HEADERS = \
 unur_distributions.h \
 unur_distributions_source.h \
 unur_stddistr.h

EXTRA_DIST = \
 README

# clean backup files
CLEANFILES = *~

# clean generated files
MAINTAINERCLEANFILES = \
 Makefile.in
