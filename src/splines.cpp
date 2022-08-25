#include <splines.h>

namespace ssim {

  // mat qr_q(const mat& X, double tol) {
  //   Rcpp::NumericMatrix nmX = Rcpp::as<Rcpp::NumericMatrix>(Rcpp::wrap(X));
  //   Rcpp::NumericMatrix nmQ = qr_q(nmX, tol);
  //   return Rcpp::as<mat>(Rcpp::wrap(nmQ));
  // }
      
  SplineBasis::SplineBasis(int order) : order(order) {
    ordm1 = order - 1;
    rdel = arma::vec(ordm1);
    ldel = arma::vec(ordm1);
    a = arma::vec(order);
  }
  SplineBasis::SplineBasis(arma::vec knots, int order) : order(order), knots(knots) {
    ordm1 = order - 1;
    nknots = knots.size();
    ncoef = nknots - order;
    rdel = arma::vec(ordm1);
    ldel = arma::vec(ordm1);
    a = arma::vec(order);
  }
  int SplineBasis::set_cursor(double x)
  {
    int i;
    /* don't assume x's are sorted */
    curs = -1; /* Wall */
    boundary = 0;
    for (i = 0; i < nknots; i++) {
      if (knots(i) >= x) curs = i;
      if (knots(i) > x) break;
    }
    if (curs > nknots - order) {
      int lastLegit = nknots - order;
      if (x == knots(lastLegit)) {
	boundary = 1; curs = lastLegit;
      }
    }
    return curs;
  }
  void
  SplineBasis::diff_table(double x, int ndiff)
  {
    int i;
    for (i = 0; i < ndiff; i++) {
      rdel(i) = knots(curs + i) - x;
      ldel(i) = x - knots(curs - (i + 1));
    }
  }
  double SplineBasis::slow_evaluate(double x, int nder)
  {
    int ti = curs, 
      lpt, apt, rpt, inner,
      outer = ordm1;
    if (boundary && nder == ordm1) { /* value is arbitrary */
      return double(0);
    }
    while(nder--) {  // FIXME: divides by zero
      for(inner = outer, apt = 0, lpt = ti - outer; inner--; apt++, lpt++)
	a(apt) = double(outer) * (a(apt + 1) - a(apt))/(knots(lpt + outer) - knots(lpt));
      outer--;
    }
    diff_table(x, outer);
    while(outer--)
      for(apt = 0, lpt = outer, rpt = 0, inner = outer + 1;
	  inner--; lpt--, rpt++, apt++)
	// FIXME: divides by zero
	a(apt) = (a(apt + 1) * ldel(lpt) + a(apt) * rdel(rpt))/(rdel(rpt) + ldel(lpt));
    return a(0);
  }
  /* fast evaluation of basis functions */
  arma::vec SplineBasis::basis_funcs(double x)
  {
    arma::vec b(order);
    diff_table(x, ordm1);
    b(0) = double(1);
    for (size_t j = 1; j <= (size_t)ordm1; j++) {
      double saved = double(0);
      for (size_t r = 0; r < j; r++) { // do not divide by zero
	double den = rdel(r) + ldel(j - 1 - r);
	if(den != double(0)) {
	  double term = b(r)/den;
	  b(r) = saved + rdel(r) * term;
	  saved = ldel(j - 1 - r) * term;
	} else {
	  if(r != double(0) || rdel(r) != double(0))
	    b(r) = saved;
	  saved = double(0);
	}
      }
      b(j) = saved;
    }
    return b;
  }
  arma:: vec SplineBasis::eval(double x, int ders) {
    arma::vec val(ncoef);
    val = arma::zeros(ncoef);
    set_cursor(x);
    int io = curs - order;
    if (io < 0 || io > nknots) {
      for (size_t j = 0; j < (size_t)order; j++) {
	val(j+io) = double(0); // R_NaN;
      }
    } else if (ders > 0) { /* slow method for derivatives */
      for(size_t i = 0; i < (size_t)order; i++) {
	for(size_t j = 0; j < (size_t)order; j++) a(j) = double(0);
	a(i) = double(1);
	val(i+io) = slow_evaluate(x, ders);
      }
    } else { 		/* fast method for value */
      arma::vec valtmp = basis_funcs(x);
      for (size_t i=0; i<valtmp.size(); i++)
	val(i+io)=valtmp(i);
    }
    return val;
  }
  arma::mat SplineBasis::basis(arma::vec x, int ders) {
    arma::mat m(x.size(), ncoef);
    for (size_t i=0; i<x.size(); i++) {
      arma::vec v = eval(x(i), ders);
      for  (size_t j=0; j<v.size(); j++)
	m(i,j)=v(j);
    }
    return m;
  }
  
  bs::bs() {} // default constructor
  bs::bs(arma::vec boundary_knots, arma::vec interior_knots, int intercept) :
    SplineBasis(4), boundary_knots(boundary_knots), interior_knots(interior_knots),
    intercept(intercept) {
    df = intercept + 3 + interior_knots.size();
    this->nknots = interior_knots.size()+8;
    this->ncoef = this->nknots - this->order;
    this->knots = arma::vec(this->nknots);
    for(size_t i=0; i<4;i++) {
      this->knots(i)=boundary_knots(0);
      this->knots(this->nknots-i-1)=boundary_knots(1);
    }
    if (interior_knots.size() > 0) 
      for(size_t i=0; i<interior_knots.size();i++) 
	this->knots(i+4)=interior_knots(i);
  }      
  arma::vec bs::eval(double x, int ders) {
    arma::vec v;
    if (x<boundary_knots(0)) {
      double k_pivot = double(0.75)*boundary_knots(0)+double(0.25)*interior_knots(0);
      double delta = x - k_pivot;
      v = bs::eval(k_pivot,0) +
	bs::eval(k_pivot,1)*delta +
	bs::eval(k_pivot,2)*delta*delta/2. +
	bs::eval(k_pivot,3)*delta*delta*delta/6.;
    }
    else if (x>boundary_knots(1)) {
      double k_pivot = double(0.75)*boundary_knots(1)+double(0.25)*interior_knots(interior_knots.size()-1);
      double delta = x - k_pivot;
      v = bs::eval(k_pivot,0) +
	bs::eval(k_pivot,1)*delta +
	bs::eval(k_pivot,2)*delta*delta/2. +
	bs::eval(k_pivot,3)*delta*delta*delta/6.;
    }
    else  {
      v = SplineBasis::eval(x, ders).subvec(1-intercept,df-intercept);
    }
    return v;
  }
  arma::mat bs::basis(arma::vec x, int ders) {
    arma::mat m(x.size(), df);
    for (size_t i=0; i<x.size(); i++) {
      arma::vec v = bs::eval(x(i), ders);
      for  (size_t j=0; j<v.size(); j++)
	m(i,j)=v(j);
    }
    return m;
  }

  ns::ns() {} // default constructor
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
  ns::ns(arma::vec boundary_knots, arma::vec interior_knots, arma::mat _q_matrix,
	 int intercept, int cure) :
    bs(boundary_knots, interior_knots, intercept), q_matrix(_q_matrix), cure(cure) {
    if (q_matrix.n_cols < q_matrix.n_rows)
      q_matrix = q_matrix.t();
    tl0 = q_matrix * bs::eval(boundary_knots(0), 0);
    tl1 = q_matrix * bs::eval(boundary_knots(0), 1);
    tr0 = q_matrix * bs::eval(boundary_knots(1), 0);
    tr1 = q_matrix * bs::eval(boundary_knots(1), 1);
  }
  arma::vec ns::eval(double x, int der) {
    if(x < this->boundary_knots(0)) {
      if (der==0) 
	return tl0 + (x - this->boundary_knots(0))*tl1;
      else if (der==1)
	return tl1;
      else return tl1*double(0);
    } else if (x > this->boundary_knots(1)) {
      if (der==0)
	return tr0 + (x - this->boundary_knots(1))*tr1;
      else if (der==1)
	return tr1;
      else return tr1*double(0);
    }
    else return q_matrix * bs::eval(x,der);
  }
  arma::mat ns::basis(arma::vec x, int ders) {
    arma::mat m(x.size(), this->df-2-cure);
    for (size_t i=0; i<x.size(); i++) {
      arma::vec v = ns::eval(x(i), ders);
      for  (size_t j=0; j<v.size(); j++)
	m(i,j)=v(j);
    }
    return m;
  }

} // namespace ssim
