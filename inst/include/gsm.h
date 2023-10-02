#ifndef MICROSIMULATION_GSM_H
#define MICROSIMULATION_GSM_H

#include <RcppArmadillo.h>
#include <splines.h>

namespace ssim {

  enum link_types {PH};

  class gsm_term {
  public:
    ns ns1;
    arma::vec gamma, x;
  };
  
  class gsm {
  public:
    link_types link_type;
    double tmin, tmax, target, target0, t0;
    arma::vec etap, etap0;
    std::vector<gsm_term> terms;
    int index;
    bool log_time;
    double link(double S);
    double linkinv(double eta);
    gsm(); // default constructor
    gsm(SEXP args);
    gsm(Rcpp::List list);
    double eta(double y);
    double eta0(double y);
    double operator()(double y);
    double rand(double tentry=0.0, int index = 0, double scale=10.0);
    double randU(double u, double tentry=0.0, int index = 0, double scale=10.0);
    double randU0(double u, int index = 0, double scale=10.0);
  };

}

#endif /* gsm.h */
