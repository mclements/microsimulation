#include <microsimulation.h>
#include <float.h> /* DBL_EPSILON */

namespace ssim {

  double rweibullHR(double shape, double scale, double hr){
    return R::rweibull(shape, scale*pow(hr,1.0/shape));
  }

  Time now() {
    return Sim::clock();
  }

  Time simTime() {
    return Sim::clock();
  }


  static Rng * default_stream, * current_stream;
  static double rn = 0.0;

  Rng::~Rng() {
    if (current_stream->id == this->id)
      current_stream = default_stream;
  }

  void Rng::set() {
    current_stream = this;
  }

  extern "C" {

    void r_create_current_stream()
    {
      default_stream = new Rng();
      current_stream = default_stream;
    }

    void r_remove_current_stream()
    {
      delete default_stream;
    }

    void r_set_user_random_seed(double * inseed) {
      double seed[6];
      for(int i=0; i<6; i++) {
	seed[i] = inseed[i];
      }
      Rng::SetPackageSeed(seed);
      default_stream->SetSeed(seed);
    }

    void r_get_user_random_seed(double * outseed) {
      double seed[6];
      default_stream->GetState(seed);
      for(int i=0; i<6; i++) {
	outseed[i] = (double)seed[i];
      }
    }

    void r_next_rng_substream() {
      default_stream->ResetNextSubstream();
    }

    void r_rng_advance_substream(double * inoutseed, int * n) {
      RngStream r;
      double seed[6];
      for (int i=0; i<6; i++)
	seed[i]=inoutseed[i];
      r.SetSeed(seed);
      r.AdvanceSubstream(0, *n);
      r.GetState(seed);
      for (int i=0; i<6; i++)
	inoutseed[i]= seed[i];
    }

    double *user_unif_rand ()
    {
      if (!current_stream) {
	REprintf("user_unif_rand(): No stream created yet!");
	return NULL;
      }
      rn = current_stream->RandU01();
      return &rn;
    }

    void test_rstream2(double * x) {
      Rng * s1 = new Rng();
      Rng * s2 = new Rng();
      x[0]=WithRNG(s1,R::rexp(1.0));
      x[1]=WithRNG(s2,R::rexp(1.0));
      s1->ResetNextSubstream();
      x[2]=R::rexp(1.0);
      delete s1;
      delete s2;
    }

  } // extern "C"

} // namespace ssim

namespace R {
  double rnormPos(double mean, double sd) {
    double x;
    while ((x=R::rnorm(mean,sd))<0.0) { }
    return x;
  }

  double rllogis(double shape, double scale) {
    double u = R::runif(0.0,1.0);
    return scale*exp(-log(1.0/u-1.0)/shape);
  }

  double rllogis_trunc(double shape, double scale, double left) {
    double S0 = 1.0/(1.0+exp(log(left/scale)*shape));
    double u = R::runif(0.0,1.0);
    return scale*exp(log(1.0/(u*S0)-1.0)/shape);
  }

  double rgompertz(double shape, double rate) {
    double u = 1.0 - R::runif(0.0, 1.0);
    return (shape < 0.0 && u<exp(rate/shape)) ? R_PosInf :
      log(1.0 - shape*log(u)/rate)/shape;
  }
}


namespace ssim {

  double R_zeroin2(			/* An estimate of the root */
					double ax,				/* Left border | of the range	*/
					double bx,				/* Right border| the root is seeked*/
					double fa, double fb,		/* f(a), f(b) */
					double (*f)(double x, void *info),	/* Function under investigation	*/
					void *info,				/* Add'l info passed on to f	*/
					double *Tol,			/* Acceptable tolerance		*/
					int *Maxit)				/* Max # of iterations */
  {
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
      *Tol = 0.0;
      *Maxit = 0;
      return a;
    }
    if(fb ==  0.0) {
      *Tol = 0.0;
      *Maxit = 0;
      return b;
    }

    while(maxit--)		/* Main iteration loop	*/
      {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	  {				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	  }
	tol_act = 2*DBL_EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	  {
	    *Maxit -= maxit;
	    *Tol = fabs(c-b);
	    return b;			/* Acceptable approx. is found	*/
	  }

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	  double t1,cb,t2;
	  cb = c-b;
	  if( a==c ) {		/* If we have only two distinct	*/
	    /* points linear interpolation	*/
	    t1 = fb/fa;		/* can only be applied		*/
	    p = cb*t1;
	    q = 1.0 - t1;
	  }
	  else {			/* Quadric inverse interpolation*/

	    q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
	    p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
	    q = (q-1.0) * (t1-1.0) * (t2-1.0);
	  }
	  if( p>(double)0 )		/* p was calculated with the */
	    q = -q;			/* opposite sign; make p positive */
	  else			/* and assign possible minus to	*/
	    p = -p;			/* q				*/

	  if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
	      && p < fabs(prev_step*q/2) )	/* and isn't too large	*/
	    new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	  if( new_step > (double)0 )	/* than tolerance		*/
	    new_step = tol_act;
	  else
	    new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	  /* Adjust c for it to have a sign opposite to that of b */
	  c = a;  fc = fa;
	}

      }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
  }

  extern "C" {
  
    void R_qMVK(int *n, double *p, double *A, double *B, double *delta, double *upper,
		double *out) {
      for (int i=0; i<*n; i++) {
	out[i] = qMVK(p[i], A[i], B[i], delta[i], upper[i]);
      }
    } 
    void R_rMVK(int *n, double *A, double *B, double *delta, double *upper,
		double *out) {
      GetRNGstate();
      for (int i=0; i<*n; i++) {
	double u = R::runif(0.0,1.0);
	out[i] = qMVK(u, *A, *B, *delta, *upper);
      }
      PutRNGstate();
    } 
    void R_rMVK2(int *n, double *A, double *B, double *delta, double *x0, double *out) {
      GetRNGstate();
      for (int i=0; i<*n; i++) {
	double u = R::runif(0.0,1.0);
	out[i] = qMVK_fast_unstable(u, *A, *B, *delta, *x0);
      }
      PutRNGstate();
    }
  }
    
}
