/*******************************************************************\
 *                                                                 *
 *   UNU.RAN -- Universal Non-Uniform Random number generator      *
 *                                                                 *
 *******************************************************************
 *   Copyright (c) 2000-2011 Wolfgang Hoermann and Josef Leydold   *
 *   Department of Statistics and Mathematics, WU Wien, Austria    *
\*******************************************************************/
#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
#  define __BEGIN_DECLS extern "C" {
#  define __END_DECLS }
#else
#  define __BEGIN_DECLS /* empty */
#  define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

#ifndef UNURAN_URNG_RNGSTREAMS_H_SEEN
#define UNURAN_URNG_RNGSTREAMS_H_SEEN

#include <RngStream.h>
UNUR_URNG *unur_urng_rngstream_new( const char *urngstr );
UNUR_URNG *unur_urng_rngstreamptr_new( RngStream rngstream );
#endif  /* UNURAN_URNG_RNGSTREAMS_H_SEEN */
__END_DECLS
