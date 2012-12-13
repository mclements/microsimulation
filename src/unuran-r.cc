     /* ------------------------------------------------------------- */
     /* File: example_rngstreams.c                                    */
     /* ------------------------------------------------------------- */

     /* ------------------------------------------------------------- */
     /* This example makes use of the RNGSTREAM library for           */
     /* for generating uniform random numbers.                        */
     /* (see http://statmath.wu.ac.at/software/RngStreams/)           */
     /* To compile this example you must have set                     */
     /*   ./configure --with-urng-rngstream                           */
     /* (Of course the executable has to be linked against the        */
     /* RNGSTREAM library.)                                           */

     /* gcc -o unuran-2 unuran-2.c -lunuran -lrngstreams -lm          */

     /* ------------------------------------------------------------- */
     
     /* Include UNURAN header files.                                  */
     #include <unuran-patched/unuran.h>
     #include <unuran-patched/uniform/unuran_urng_rngstreams.h>

     #include <string>
     #include <cstdarg>
     
     /* ------------------------------------------------------------- */

using std::string;

// http://stackoverflow.com/questions/2342162/stdstring-formatting-like-sprintf
std::string string_format(const std::string &fmt, ...) {
    int size=100;
    std::string str;
    va_list ap;
    while (1) {
        str.resize(size);
        va_start(ap, fmt);
        int n = vsnprintf((char *)str.c_str(), size, fmt.c_str(), ap);
        va_end(ap);
        if (n > -1 && n < size) {
            str.resize(n);
            return str;
        }
        if (n > -1)
            size=n+1;
        else
            size*=2;
    }
}

class UnuranGen {
public:
  UNUR_GEN * g;
  UnuranGen(UNUR_GEN * gin = 0) : g(gin) {}; // default constructor?
  virtual ~UnuranGen() { 
    if (g != NULL)
      unur_free(g); 
  }
  virtual double sample() = 0;
};

class UnuranCont : public UnuranGen {
public:
  virtual double sample() { 
    return unur_sample_cont(g);
  }
};
  
class UnuranDiscr : public UnuranGen {
public:
  virtual double sample() { 
    return unur_sample_discr(g);
  }
};

class Unuran {
public:
  UNUR_URNG *urng;
  Unuran(const char * urngstr) {
    urng = unur_urng_rngstream_new(urngstr);
    if (urng == NULL) exit (EXIT_FAILURE);
  }
  ~Unuran() {
    if (urng != NULL)
      unur_urng_free(urng);
  }
  void free() {
    if (urng != NULL)
      unur_urng_free(urng);
  }
  void nextsub() {
    unur_urng_nextsub(urng);
  }
  UNUR_GEN *gen(const string s) {
    UNUR_GEN * g = unur_str2gen(s.c_str());
    if (g == NULL) {
      fprintf(stderr, "ERROR: cannot create generator object\n");
      exit (EXIT_FAILURE);
    }    
    unur_chg_urng(g, urng);
    return g;
  }
  UNUR_GEN * urbeta(double shape1, double shape2, double lb = 0.0, double ub = 1.0) {
    string s = string_format("beta(%f,%f); domain=(%f,%f) & method=HINV",
			     shape1, shape2, lb, ub);
    return gen(s);
  }
  UNUR_GEN *urbinom(int size, double prob, int lb = 0, int ub = -1) {
    string s = string_format("binomial(%d,%f); domain=(%d,%d)", 
			     size, prob, lb, ub == -1 ? size : ub);
    UNUR_DISTR* distr = unur_str2distr(s.c_str());
    unur_distr_discr_make_pv(distr); 
    UNUR_PAR* par = unur_dgt_new (distr);
    UNUR_GEN* gen = unur_init(par); // this destroys the par object
    if (gen == NULL) {
      fprintf(stderr, "ERROR: cannot create generator object\n");
      exit (EXIT_FAILURE);
    }    
    unur_chg_urng(gen, urng);
    unur_distr_free(distr);
    return gen;
  }
  UNUR_GEN *urexp(double rate = 1.0, double lb = 0.0, double ub = UNUR_INFINITY) {
    string s = string_format("exponential(%f); domain=(%f,%f) & method=CSTD",
			     1.0/rate, lb, ub);
    return gen(s);
  }
  UNUR_GEN *urgamma(double shape, double scale = 1.0, double lb = 0.0, 
		   double ub = UNUR_INFINITY) {
    string s = string_format("gamma(%f,%f); domain=(%f,%f) & method=HINV",
			     shape, scale, lb, ub);
    return gen(s);
  }
  UNUR_GEN *urlnorm(double meanlog = 0.0, double sdlog = 1.0, double lb = 0.0, 
		   double ub = UNUR_INFINITY) {
    string s = string_format("lognormal(%f,%f,%f); domain=(%f,%f) & method=HINV",
			     meanlog, sdlog, lb, lb, ub);
    return gen(s);
  }
  UNUR_GEN *urnorm(double mean = 0.0, double sd = 1.0, double lb = -UNUR_INFINITY, 
		   double ub = UNUR_INFINITY) {
    string s = string_format("normal(%f,%f); domain=(%f,%f) & method=HINV",
			     mean, sd, lb, ub);
    return gen(s);
  }
  UNUR_GEN *urweibull(double shape, double scale = 1.0, double lb = 0.0, double ub = UNUR_INFINITY) {
    string s = string_format("weibull(%f,%f); domain=(%f,%f) & method=HINV",
			     shape, scale, lb, ub);
    return gen(s);
  }
};


double sampleOnce(UNUR_GEN *g) {
  double out = unur_sample_cont(g);
  unur_free(g);
  return out;
}
     
