#include <R_ext/Applic.h>
#include <Rcpp.h>

typedef double optimfn(int, double *, void *);
typedef void optimgr(int, double *, double *, void *);

// void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fminfn,
// 	   int *fail, double abstol, double intol, void *ex,
// 	   double alpha, double bet, double gamm, int trace,
// 	   int *fncount, int maxit)

class NelderMead {
public:
  NelderMead(int trace = 0, int maxit = 500, 
	     double abstol = - INFINITY,
	     double reltol = 1.0e-8, 
	     double alpha = 1.0, double beta = 0.5, double gamma = 2.0) : 
    trace(trace), maxit(maxit), abstol(abstol), reltol(reltol), 
    alpha(alpha), beta(beta), gamma(gamma) { 
  }
  void optim(optimfn fn, Rcpp::NumericVector init, Rcpp::NumericVector X, void * ex) {
    n = init.size();
    nmmin(n, &init[0], &X[0], &Fmin, fn,
	  &fail, abstol, reltol, ex,
	  alpha, beta, gamma, trace,
	  &fncount, maxit);
  }
  int n, trace, maxit, fail, fncount;
  double abstol, reltol, alpha, beta, gamma, Fmin;
};

// void
// vmmin(int n0, double *b, double *Fmin, optimfn fminfn, optimgr fmingr,
//       int maxit, int trace, int *mask,
//       double abstol, double reltol, int nREPORT, void *ex,
//       int *fncount, int *grcount, int *fail)

class BFGS {
public:
  BFGS(int trace = 0, int maxit = 100, 
       double abstol = - INFINITY,
       double reltol = 1.0e-8, int report = 10) : 
    trace(trace), maxit(maxit), report(report), abstol(abstol), reltol(reltol) { }
  void optim(optimfn fn, optimgr gr, Rcpp::NumericVector init, void * ex) {
    n = init.size();
    std::vector<int> mask(n,1);
    vmmin(n, &init[0], &Fmin, fn, gr, maxit, trace, &mask[0], abstol, reltol, report,
	  ex, &fncount, &grcount, &fail);
  }
  int n, trace, maxit, report, fncount, grcount, fail;
  double abstol, reltol, Fmin;
};


double fminfn(int n, double * Bvec, void *ex) {
  double beta = * (double *) ex;
  return Bvec[0]*Bvec[0] + (Bvec[1]- beta)*(Bvec[1]- beta);
}

void grfn(int n, double * Bvec, double * gr, void *ex) {
  double beta = * (double *) ex;
  gr[0] = 2*Bvec[0];
  gr[1] = 2*(Bvec[1]- beta);
}


RcppExport SEXP test_nmmin(SEXP Sinit) {

  double ex = 1.0;
  
  Rcpp::NumericVector init = Rcpp::as<Rcpp::NumericVector>(Sinit);
  Rcpp::NumericVector X(init.size());
  Rcpp::NumericVector X2(Rcpp::clone(init));

  NelderMead nm;
  BFGS bfgs;
  
  nm.optim(fminfn, init, X, (void *) &ex);

  Rprintf("X[0]=%g\nX[1]=%g\n",X[0],X[1]);
  Rprintf("fail=%i\nFmin=%g\nfncount=%i\n",nm.fail,nm.Fmin,nm.fncount);

  bfgs.optim(fminfn, grfn, X2, (void *) &ex);

  Rprintf("X[0]=%g\nX[1]=%g\n",X2[0],X2[1]);
  Rprintf("fail=%i\nFmin=%g\nfncount=%i\ngrcount=%i\n",bfgs.fail,bfgs.Fmin,bfgs.fncount,bfgs.grcount);

  return Rcpp::wrap(X2);

}
// R CMD INSTALL ~/src/R/microsimulation
// R -q -e "require(microsimulation); .Call('test_nmmin',1:2,PACKAGE='microsimulation')"
