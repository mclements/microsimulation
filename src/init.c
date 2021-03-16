#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void r_get_user_random_seed(void *);
extern void r_next_rng_substream();
extern void r_rng_advance_substream(void *, void *);
extern void r_set_user_random_seed(void *);
extern void r_create_current_stream();
extern void r_remove_current_stream();

/* .Call calls */
extern SEXP callCalibrationSimulation(SEXP);
extern SEXP callIllnessDeath(SEXP);
extern SEXP callPersonSimulation(SEXP, SEXP);
extern SEXP callSimplePerson(SEXP);
extern SEXP callSimplePerson2(SEXP);
extern SEXP pqueue__cancel(SEXP, SEXP);
extern SEXP pqueue__clear(SEXP);
extern SEXP pqueue__empty(SEXP);
extern SEXP pqueue__new(SEXP);
extern SEXP pqueue__pop(SEXP);
extern SEXP pqueue__push(SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"r_get_user_random_seed",  (DL_FUNC) &r_get_user_random_seed,  1},
    {"r_next_rng_substream",    (DL_FUNC) &r_next_rng_substream,    0},
    {"r_rng_advance_substream", (DL_FUNC) &r_rng_advance_substream, 2},
    {"r_set_user_random_seed",  (DL_FUNC) &r_set_user_random_seed,  1},
    {"r_create_current_stream",   (DL_FUNC) &r_create_current_stream,   0},
    {"r_remove_current_stream",   (DL_FUNC) &r_remove_current_stream,   0},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"callCalibrationSimulation", (DL_FUNC) &callCalibrationSimulation, 1},
    {".callIllnessDeath",          (DL_FUNC) &callIllnessDeath,          1},
    {".callPersonSimulation",      (DL_FUNC) &callPersonSimulation,      2},
    {".callSimplePerson",          (DL_FUNC) &callSimplePerson,          1},
    {".callSimplePerson2",         (DL_FUNC) &callSimplePerson2,         1},
    {"pqueue__cancel",            (DL_FUNC) &pqueue__cancel,            2},
    {"pqueue__clear",             (DL_FUNC) &pqueue__clear,             1},
    {"pqueue__empty",             (DL_FUNC) &pqueue__empty,             1},
    {"pqueue__new",               (DL_FUNC) &pqueue__new,               1},
    {"pqueue__pop",               (DL_FUNC) &pqueue__pop,               1},
    {"pqueue__push",              (DL_FUNC) &pqueue__push,              3},
    {NULL, NULL, 0}
};

void R_init_microsimulation(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
