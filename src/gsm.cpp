#include <RcppArmadillo.h>
#include <splines.h>
#include <c_optim.h>
#include <gsm.h>

namespace ssim {

  double gsm::link(double S) {
    return link_type==PH ? std::log(-std::log(S)) : -100.0;
  }
  double gsm::linkinv(double eta) {
    return link_type==PH ? std::exp(-std::exp(eta)) : 1.0e-10;
  }
  gsm::gsm() {}
  double gsm::eta(double y) {
    double eta = etap(index);
    for (size_t i=0; i<terms.size(); i++)
      if (terms[i].x(index) != 0.0)
	eta += terms[i].x(index) * arma::sum(terms[i].ns1.eval(y,0) % terms[i].gamma);
    return eta;
  }
  double gsm::operator()(double y) {
    return eta(y) - target;
  }
  double gsm::rand(double tentry, int index, double scale) {
    Rcpp::RNGScope rngScope;
    using std::log;
    double u = R::runif(0.0,1.0);
    return randU(u, tentry, index, scale);
  }
  double gsm::randU(double u, double tentry, int index, double scale) {
    using std::log;
    double ymin = tentry == 0.0 ? (log_time ? log(tmin/scale) : tmin/scale) : (log_time ? log(tentry) : tentry);
    double ymax = log_time ? log(tmax*scale) : tmax*scale;
    this->index = index;
    target = (tentry==0.0 ? link(u) : link(u*linkinv(eta(ymin))));
    double root = std::get<0>(R_zeroin2_functor_ptr<gsm>(ymin, ymax, this, 1.0e-8, 100));
    return log_time ? std::exp(root) : root;
  }

  gsm::gsm(Rcpp::List list) {
    try {
      using namespace Rcpp;
      std::string link_name = as<std::string>(list("link_name"));
      tmin = as<double>(list("tmin"));
      tmax = as<double>(list("tmax"));
      double inflate = as<double>(list("inflate"));
      tmin = tmin/inflate; tmax = tmax*inflate;
      etap = as<arma::vec>(list("etap"));
      List lterms = as<List>(list("terms"));
      for (int i=0; i<lterms.size(); i++) {
	List lterm = as<List>(lterms(i));
	gsm_term term;
	term.gamma = as<arma::vec>(lterm("gamma"));
	arma::vec knots = as<arma::vec>(lterm("knots"));
	arma::vec Boundary_knots = as<arma::vec>(lterm("Boundary_knots"));
	int intercept = as<int>(lterm("intercept"));
	arma::mat q_const = as<arma::mat>(lterm("q_const"));
	int cure = as<int>(lterm("cure"));
	term.ns1 = ns(Boundary_knots, knots, q_const, intercept, cure);
	term.x = as<arma::vec>(lterm("x"));
	terms.push_back(term);
      }
      log_time = as<bool>(list("log_time"));
      target = 0.0;
      index = 0;
      if (link_name == "PH") link_type = PH;
    } catch(std::exception &ex) {	
      forward_exception_to_r(ex);
    } catch(...) { 
      ::Rf_error("c++ exception (unknown reason)"); 
    }
  } 

  gsm::gsm(SEXP args) : gsm(Rcpp::as<Rcpp::List>(args)) { }
  
  RcppExport SEXP test_read_gsm(SEXP gsm_args, SEXP start_args) {
    Rcpp::RNGScope rngScope;
    double start = Rcpp::as<double>(start_args);
    gsm gsm1(gsm_args);
    return Rcpp::wrap(gsm1.rand(start));
  }
  
} // namespace ssim
