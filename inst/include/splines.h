#ifndef SPLINES_H
#define SPLINES_H

#include <RcppArmadillo.h>

namespace ssim {

  /* arma::mat qr_q(const arma::mat& X, double tol = 1E-12); */
      
  class SplineBasis {
  public:
    int order,			/* order of the spline */
      ordm1,			/* order - 1 (3 for cubic splines) */
      nknots,			/* number of knots */
      curs,			/* current position in knots vector */
      boundary,		/* must have knots[(curs) <= x < knots(curs+1) */
      ncoef;			/* number of coefficients */
    /* except for the boundary case */
    arma::vec ldel;  	/* differences from knots on the left */
    arma::vec rdel;	/* differences from knots on the right */
    arma::vec knots;	/* knot vector */
    arma::vec coeff;	/* coefficients */
    arma::vec a;		/* scratch array */
    SplineBasis(int order = 4);
    SplineBasis(arma::vec knots, int order = 4);
    int set_cursor(double x);
    void diff_table(double x, int ndiff);
    double slow_evaluate(double x, int nder);
    /* fast evaluation of basis functions */
    arma::vec basis_funcs(double x);
    arma::vec eval(double x, int ders=0);
    arma::mat basis(arma::vec x, int ders=0);
  };
  
  class bs : public SplineBasis {
  public:
    arma::vec boundary_knots, interior_knots;
    int intercept, df;
    bs(); // default constructor
    bs(arma::vec boundary_knots, arma::vec interior_knots, int intercept = 0);
    arma::vec eval(double x, int ders=0);
    arma::mat basis(arma::vec x, int ders=0);
  };

  class ns : public bs {
  public:
    arma::vec tl0, tl1, tr0, tr1;
    arma::mat q_matrix;
    int cure;
    ns(); // default constructor
    // ns(vec boundary_knots, vec interior_knots, int intercept=0) :
    //   bs(boundary_knots, interior_knots, intercept) {
    //   // calculate the Q matrix
    //   mat const_basis = bs::basis(boundary_knots, 2);
    //   mat qd = qr_q(const_basis.t());
    //   mat qsub(qd.n_rows, qd.n_cols-2);
    //   for (size_t i=0; i<qsub.n_rows; i++)
    // 	for (size_t j=0; j<qsub.n_cols; j++)
    // 	  qsub(i,j) = qd(i,j+2);
    //   q_matrix = qsub.t();
    //   tl0 = q_matrix * bs::eval(boundary_knots(0), 0);
    //   tl1 = q_matrix * bs::eval(boundary_knots(0), 1);
    //   tr0 = q_matrix * bs::eval(boundary_knots(1), 0);
    //   tr1 = q_matrix * bs::eval(boundary_knots(1), 1);
    // }
    ns(arma::vec boundary_knots, arma::vec interior_knots, arma::mat _q_matrix,
       int intercept=0, int cure=0);
    arma::vec eval(double x, int der);
    arma::mat basis(arma::vec x, int ders=0);
  }; // class ns

} // namespace ssim

#endif /* splines.h */
