refresh
require(foreign)
require(rstpm2)
require(dplyr)
require(ggplot2)
require(lattice)
require(maxLik)
##
temp <- read.dta("~/Downloads/rescreen-20141024.dta")
temp2 <- temp %>%
    filter((is.na(diadate) | sample_date<diadate) & time>0) %>% # restrict to non-cancers to PSA before cancer
    mutate(total_cat=cut(total,c(0,1,3,10,Inf))) # generate PSA categories
temp4 <- temp2 %>%
    mutate(age_start = as.numeric((sample_date-dob)/365),
           age_finish = as.numeric((pmin(next_sample_date,finish,na.rm=TRUE)-dob)/365)) %>%
    as.data.frame()
temp4 <- filter(temp4, !is.na(total_cat)) %>%
    mutate(time = time/365.25)
temp2b <- read.dta("/home/marcle/Downloads/psa_rates_nopca_20141017.dta") # (#test | no pca)
## xyplot(psa_rate ~ year | factor(age5), data=temp2b, subset=age5>=40 & year<=2012, type="l")
## Results are moderately stable for 2008-2011
expit <- rstpm2:::expit
Call <- function(args,fun) do.call(deparse(substitute(fun)),as.list(args))
##

cureModel <- function(formula,data=NULL,formula.shape=~1,formula.cure=~1,par=NULL,maxLik=FALSE,...) {
    formula.scale <- formula
    rstpm2:::lhs(formula.scale) <- NULL
    mf <- model.frame(formula,data)
    y <- model.extract(mf,"response")
    if (ncol(y)==2) { # right censored
        delayed <- FALSE
        time <- y[,1]
        status <- y[,2]
    }
    if (ncol(y)==3) { # left truncated and right censored
        entry <- y[,1]
        time <- y[,2]
        status <- y[,3]
        if (delayed <- any(entry>0))
            i <- which(entry>0)
    }
    formulae <- list(shape=formula.shape,
                     scale=formula.scale,
                     cure=formula.cure)
    X <- lapply(formulae,model.matrix,data)
    stopifnot(nrow(X$scale) == nrow(X$shape))
    n <- lapply(X,ncol)
    cumn <- cumsum(n)
    index <- mapply(function(a,b) (a+1):b, c(0,cumn[-length(cumn)]),cumn,SIMPLIFY=FALSE)
    names(index) <- names(formulae)
    ll <-function(beta) {
        lp <- mapply(function(x,idx) as.vector(x %*% beta[idx]), X, index, SIMPLIFY=FALSE)
        shape <- exp(lp$shape)
        scale <- exp(lp$scale)
        cure <- rstpm2:::expit(lp$cure)
        S <- cure+(1-cure)*pweibull(time,shape,scale,lower.tail=FALSE)
        f <- (1-cure)*dweibull(time,shape,scale)
        value <- sum(ifelse(status==1, log(f), log(S)))
        if (delayed) {
            S0 <- cure[i]+(1-cure[i])*pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
            value <- value - sum(log(S0))
        }
        value
    }
    ## logLik <-function(beta) {
    ##     lp <- mapply(function(x,idx) as.vector(x %*% beta[idx]), X, index, SIMPLIFY=FALSE)
    ##     shape <- exp(lp$shape)
    ##     scale <- exp(lp$scale)
    ##     cure <- rstpm2:::expit(lp$cure)
    ##     S <- cure+(1-cure)*pweibull(time,shape,scale,lower.tail=FALSE)
    ##     f <- (1-cure)*dweibull(time,shape,scale)
    ##     value <- ifelse(status==1, log(f), log(S))
    ##     if (delayed) {
    ##         S0 <- cure[i]+(1-cure[i])*pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
    ##         value[i] <- value[i] + log(S0)
    ##     }
    ##     value
    ## }
    if (is.null(par)) {
        aft <- survreg(formula,data)
        ## AFT -> Weibull
        sigma <- aft$scale # AFT "scale"
        shape <- 1/sigma # Weibull "shape"
        logscale <- coef(aft)[1]
        logbeta <- (coef(aft)[-1])
        par <- c(shape=c(log(shape),rep(0,n$shape-1)),
                 scale=c(logscale,logbeta),
                 cure=c(-7,rep(0,n$cure-1)))
        print(par)
        ## par <- c(shape=rep(0,n.shape), # shape
        ##                scale=c(log(sum(status)/sum(time)),rep(0,n.scale-1)), # scale
        ##                cure=c(logit((1-mean(status))/2), rep(0,n.cure-1)))
    }
    if (maxLik) obj <- maxLik(ll, start = par, ...) else {
        fit <- .Call("fitCureModel",time,status,X$shape,X$scale,X$cure,par,package="rstpm2") # fast
        coef <- fit$coef
        obj <- list(maximum=-fit$Fmin,
                    estimate=coef,
                    gradient=rstpm2:::grad(ll,coef),
                    gradientObs=matrix(NA,length(time),length(coef)), # hack
                    hessian=fit$hessian,
                    code=if (fit$fail==0) 1 else 4,
                    message = "",
                    last.step=NULL,
                    fixed = rep(FALSE,length(coef)),
                    iterations=0,
                    type="",
                    constraints=NULL,
                    varcovar=solve(fit$hessian))
        class(obj) <- c("maxLik","maxim","list")
        ## ## alternative approach using bbmle::mle2
        ## negll <- function(beta) -ll(beta)
        ## rownames(fit$hessian) <- colnames(fit$hessian) <- parnames(negll) <- names(coef)
        ## obj <- mle2(negll, start = par, eval.only = TRUE, vec.par = TRUE) # SLOW
        ## obj@details$convergence <- fit$fail
        ## obj@vcov <- solve(fit$hessian)
    }
    obj
}
temp50 <- filter(temp4, age5==50)
## debug(cureModel)

system.time(fit <- cureModel(Surv(time,event)~1, data=temp50, par=c(shape=0.127633659438804, scale=0.730222998678452, cure=-2.0680040438727)))
system.time(fit <- cureModel(Surv(time,event)~1, data=temp50))

summary(fit2 <- cureModel(Surv(time,event)~total_cat, data=temp50))
## Log-Likelihood: -113312
summary(fit2 <- cureModel(Surv(time,event)~total_cat, data=temp50, par=c(shape=0.147551370026859, scale=c(0.884863166386332, -0.115910931108631, -1.09037223824094, -1.52389713434596), cure=-2.18162256296465)))

## Simple models:
temp30plus <- filter(temp4, age5>=30) %>% group_by(age5, total_cat)
models <- do(temp30plus,cureModel(Surv(time,event)~1, data=.) %>% coef %>% Call(data.frame))
models2 <- mutate(models, cure=expit(cure), shape=exp(shape), scale=exp(scale..Intercept.),
                  scale..Intercept.=NULL, psa=switch(unclass(total_cat),"1"=0,"2"=1,"3"=3,"4"=10)) %>%
    ungroup() %>% mutate(total_cat=psa,psa=NULL)
dput(as.data.frame(models2))

## Weibull regression formulated in terms of the *weibull functions (cf. survreg)
WeibullRegression <- function(formula,data=NULL,formula.shape=~1,par=NULL,...) {
    formula.scale <- formula
    rstpm2:::lhs(formula.scale) <- NULL
    mf <- model.frame(formula,data)
    y <- model.extract(mf,"response")
    if (ncol(y)==2) { # right censored
        delayed <- FALSE
        time <- y[,1]
        status <- y[,2]
    }
    if (ncol(y)==3) { # left truncated and right censored
        entry <- y[,1]
        time <- y[,2]
        status <- y[,3]
        if (delayed <- any(entry>0))
            i <- which(entry>0)
    }
    formulae <- list(shape=formula.shape,
                     scale=formula.scale)
    X <- lapply(formulae,model.matrix,data)
    stopifnot(nrow(X$scale) == nrow(X$shape))
    cumn <- cumsum(sapply(X,ncol))
    index <- mapply(function(a,b) (a+1):b, c(0,cumn[-length(cumn)]),cumn,SIMPLIFY=FALSE)
    names(index) <- names(formulae)
    ll <-function(beta) {
        lp <- mapply(function(x,idx) as.vector(x %*% beta[idx]), X, index, SIMPLIFY=FALSE)
        shape <- exp(lp$shape)
        scale <- exp(lp$scale)
        S <- pweibull(time,shape,scale,lower.tail=FALSE)
        f <- dweibull(time,shape,scale)
        value <- sum(ifelse(status==1, log(f), log(S)))
        if (delayed) {
            S0 <- pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
            value <- value - sum(log(S0))
        }
        value
    }
    if (is.null(par))
        par <- c(shape=rep(0,n.shape), # shape
                       scale=c(log(sum(status)/sum(time)),rep(0,n.scale-1)))
    maxLik(ll, start = par, ...)
}
system.time(fit2 <- WeibullRegression(Surv(time,event)~total_cat, data=temp50, par=c(shape=0.147551370026859, scale=c(0.884863166386332, -0.115910931108631, -1.09037223824094, -1.52389713434596))))
summary(fit2)
summary(fit.init2 <- survreg(Surv(time,event)~total_cat,data=temp50))

## AFT -> Weibull
sigma <- fit.init2$scale # AFT "scale"
shape <- 1/sigma # Weibull "shape"
logscale <- coef(fit.init2)[1]
logbeta <- (coef(fit.init2)[-1])
c(log(shape),logscale,logbeta)
summary(fit2)

## maxLik output
solve(vcov(fit))
exp(coef(fit)[1:2])
expit(coef(fit)[3])

## > solve(vcov(fit$maxLik))
##           shape     scale      cure
## shape 81839.971  4831.236 -4540.198
## scale  4831.236 56432.327  7144.990
## cure  -4540.198  7144.990  4205.504


## plot of the distributions
x <- seq(0,20,length=101)
rescreening2 <- lapply(1:nrow(rescreening),
                       function(i) data.frame(age5=rescreening$age5[i],
                                              total_cat=rescreening$total_cat[i],
                                              x=x,
                                              p=(1-rescreening$cure[i])*pweibull(x,rescreening$shape[i],rescreening$scale[i])))
rescreening2 <- do.call("rbind",rescreening2)
xyplot(p ~ x | factor(total_cat)+factor(age5),data=rescreening2,type="l",subset=age5>=55) 


fit2 <- cureModel(Surv(time,event)~total_cat, data=temp50, par=c(shape=0.147551370026859, scale=c(0.884863166386332, -0.115910931108631, -1.09037223824094, -1.52389713434596), cure=-2.18162256296465))
fit3 <- cureModel(Surv(time,event)~total_cat,
                  formula.shape=~total_cat,
                  data=temp50,
                  par=c(shape=c(0.27,-0.1,-0.4,-0.54),
                      scale=c(0.85,-0.1,-1.2,-1.8),
                      cure=-2.05)) # best fit?
fit4 <- cureModel(Surv(time,event)~total_cat,
                  formula.cure=~total_cat,
                  data=temp50,
                  par=c(0.154696663520933, 0.865091437590591, -0.111246732543863, -1.05240102943909,
-1.48716715943421, -2.04615245632196, -0.0207990798623829, -0.639229272169032,
-1.17752470254301))
fit5 <- cureModel(Surv(time,event)~total_cat,
                  formula.shape=~total_cat,
                  formula.cure=~total_cat,
                  data=temp50,
                  par=c(shape=c(0.27,-0.1,-0.4,-0.54),
                      scale=c(0.85,-0.1,-1.2,-1.8),
                      cure=c(-2.05,0,0,0)))

sapply(list(fit,fit2,fit3,fit4,fit5),AIC)


## check the screening results from the simulation
require(microsimulation)
require(sqldf)
require(mgcv)
set.seed(12345)
sim <- callFhcrc(1e6, screen="screenUptake")
events <- sim$summary$events
pt <- sim$summary$pt
rates <- sqldf('with pop as (select age, year, sum(pt) as pt from pt where year between 1980 and 2020 and age between 40 and 89 group by age, year), screen as (select age, year, count(*) as n from events where event in ("toScreen","toBiopsyFollowUpScreen") group by age, year) select pop.age, pop.year, pop.pt, coalesce(n,0) as n, 1.0*coalesce(n,0)/pt as rate from pop left natural join screen')
rates$pred <- predict(gam(n~s(age,year)+offset(log(pt)),data=rates,family="poisson"),newdata=transform(rates,pt=1),type="response")
ratestab <- xtabs(rate~age+year,rates)
## ratestab[ratestab==0] <- NA
filled.contour(as.numeric(dimnames(ratestab)$age),
               as.numeric(dimnames(ratestab)$year),
               ratestab)
ratestab <- xtabs(n~age+year,rates)
filled.contour(as.numeric(dimnames(ratestab)$age),
               as.numeric(dimnames(ratestab)$year),
               ratestab)

## Screening uptake
require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
using namespace Rcpp;
  double rllogis(double shape, double scale) {
    double u = R::runif(0.0,1.0);
    return scale*exp(-log(1.0/u-1.0)/shape);
  }
  double rllogis_trunc(double shape, double scale, double left) {
    double S0 = 1.0/(1.0+exp(log(left/scale)*shape));
    double u = R::runif(0.0,1.0);
    return scale*exp(log(1.0/(u*S0)-1.0)/shape);
  }
// [[Rcpp::export]]
NumericVector testRcpp(NumericVector cohort) {
      //double pscreening = cohort[0]>=1932.0 ? 0.9 : 0.9-(1932.0 - cohort[0])*0.01;
      double pscreening = 0.9;
      double shapeA = 3.8;
      double scaleA = 15.0;
      double shapeT = 2.16;
      double scaleT = 11.7;
      double uscreening = R::runif(0.0,1.0);
      double first_screen;
      if (cohort[0] > 1960.0) {
	first_screen = 35.0 + rllogis(shapeA,scaleA); // (i) age
      } else if (cohort[0] < 1945.0) {
	first_screen = (1995.0 - cohort[0]) + rllogis(shapeT,scaleT); // (ii) period
      } else {
	double age0 = 1995.0 - cohort[0];
	double u = R::runif(0.0,1.0);
	if ((age0 - 35.0)/15.0 < u) // (iii) mixture
	  first_screen = age0 + rllogis_trunc(shapeA,scaleA,age0-35.0);
	else first_screen = age0 + rllogis(shapeT,scaleT);
      }
      if (uscreening<pscreening)
	return wrap(first_screen);
      else return wrap(-1.0);
}
'
sourceCpp(code=src)
tmp <- sapply(rep(1950,1000),testRcpp)
mean(tmp==-1) # ~10%
i <- tmp != -1
plot(density(tmp[i]),xlim=c(30,120))

## Re-screening
require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
#include "/home/marcle/src/R/microsimulation/src/rcpp_table.h"
using namespace Rcpp;
  template<class T>
  T bounds(T x, T a, T b) {
    return (x<a)?a:((x>b)?b:x);
  }
// [[Rcpp::export]]
DataFrame testRcpp(DataFrame rescreening, double time, double psa) {
  typedef Table<pair<double,double>,double> TableDDD; // as per TableBiopsyCompliance
  TableDDD rescreen_shape, rescreen_scale, rescreen_cure;
  rescreen_shape = TableDDD(rescreening, "age5", "total", "shape");
  rescreen_scale = TableDDD(rescreening, "age5", "total", "scale");
  rescreen_cure  = TableDDD(rescreening, "age5", "total", "cure");
  TableDDD::key_type key = TableDDD::key_type(bounds<double>(time,30.0,90.0),psa);
  double prescreened = 1.0 - rescreen_cure(key);
  double shape = rescreen_shape(key);
  double scale = rescreen_scale(key);
  double u = R::runif(0.0,1.0);
  double t = time + R::rweibull(shape,scale);
  bool brescreened = u<prescreened;
  return wrap(Rcpp::DataFrame::create(_["u"]=u,_["prescreened"]=prescreened,_["t"]=t,_["shape"]=shape,_["scale"]=scale,_["brescreened"]=brescreened));
}
'
sourceCpp(code=src)
rescreening$total <- rescreening$total_cat
tmp <- do.call("rbind",lapply(1:1e3,function(i) testRcpp(rescreening,71,0.1)))
head(tmp)
mean(tmp$brescreened)
plot(density(tmp$t[tmp$brescreened],from=71))

subset(rescreening,age5==60)
subset(rescreening,age5==30)

require(Rcpp)
require(RcppArmadillo)
require(inline)
src <- '
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector testRcpp(NumericVector inX) {
vec X = as<vec>(inX);
//X.resize(X.size()-1);
//return wrap(X==math::inf() || X==-math::inf());
return wrap(X);
}
'
sourceCpp(code=src)
testRcpp(as.double(c(-Inf,1,Inf,NA)))

require(Rcpp)
require(RcppArmadillo)
require(inline)
src <- '
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericVector testRcpp(NumericVector inX) {
rowvec X = as<rowvec>(inX);
return wrap(exp(X));
}
'
sourceCpp(code=src)
testRcpp(as.double(c(-Inf,1,Inf,NA)))


require(Rcpp)
require(RcppArmadillo)
require(inline)
src <- '
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix testRcpp(NumericMatrix inX) {
mat X = as<mat>(inX);
return wrap(expmat(X));
}
'
sourceCpp(code=src)
testRcpp(rbind(c(-0.1,0.1),c(0.2,-0.2)))



require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List testRcpp() {
List list;
list["a"].push_back(1); // FAILS
return wrap(list);
}
'
sourceCpp(code=src)
testRcpp()

require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP myeval(SEXP sexp) {
Environment env = Environment::global_env();
switch(TYPEOF(sexp)) {
 case REALSXP:
 case INTSXP:
 case LGLSXP:
  return sexp;
  break;
 case SYMSXP:
  return Rf_eval(sexp,env);
  break;
 case LANGSXP:
  //
  break;
};
return wrap(0);
}
'
sourceCpp(code=src)
x=3:4
myeval(as.name("x"))
myeval(1:10)
myeval(T)

## We can readily evaluate values and expressions using R
require(Rcpp)
require(inline)
src <- '#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
SEXP myeval(SEXP sexp) {
  Environment env = Environment::global_env();
  return Rf_eval(sexp,env);
}'
sourceCpp(code=src)
x=3:4
myeval(as.name("x"))
myeval(1:10)
myeval(T)
myeval(exp(x))



require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
// [[Rcpp::export]]
  double rnormPos1(double mean, double sd) {
    double x;
    while ((x=R::rnorm(mean,sd))<0.0) { }
    return x;
  }
'
sourceCpp(code=src)
y <- sapply(1:1e4,function(i) rnormPos1(1,2))
mean(y)
sd(y)
rnormPos <- function(n,mean=0,sd=1,lbound=0) {
    if (length(mean)<n) mean <- rep(mean,length=n)
    if (length(sd)<n) sd <- rep(sd,length=n)
    x <- rnorm(n,mean,sd)
    while(any(i <- which(x<lbound)))
        x[i] <- rnorm(length(i),mean[i],sd[i])
    x
}
mean(rnormPos(1e5,1,2))
sd(rnormPos(1e5,1,2))
trunc_mean <- function(mu,sigma) {
    alpha <- -mu/sigma
    lambda <- dnorm(alpha)/pnorm(alpha,lower.tail=FALSE)
    mu + sigma*lambda
}
trunc_var <- function(mu,sigma) {
    alpha <- -mu/sigma
    lambda <- dnorm(alpha)/pnorm(alpha,lower.tail=FALSE)
    delta <- lambda*(lambda-alpha)
    sigma^2*(1-delta)
}
trunc_mean(1,2)
sqrt(trunc_var(1,2))

require(Rcpp)
src <- "#include <Rcpp.h>
using namespace Rcpp ;
// [[Rcpp::export]]
SEXP fn_impl( Language call, Environment env){
    for (int i=0; i<1e5; ++i)
       Rf_eval( call, env ) ;
    return Rf_eval( call, env ) ;
}"
sourceCpp(code=src)
fn <- function(call, list){
    fn_impl( call, list2env(list) )
}
fn(substitute(a*2), list(a = 1)) # fast
fn(substitute(exp(x)), list(x = 1)) # fast
## system.time(fn(substitute(lm(y~x)), list(x = 1:10, y=1:10))) # slow
system.time(for (i in 1:1e5) { x = 1:10; y=1:10; x+y })
system.time(fn(substitute(x+y), list(x = 1:10, y=1:10))) # 3-4 times faster

##
dots <- function(...) {
    eval(substitute(alist(...)))
}
##
let <- function(...) {
    args <- as.list(sys.call())[-1]
    n <- length(args)
    eval(args[[n]], args[-n])
}
let(a=1,b=2,a+b)
## let(a=1,b=2,let(c=3,a+b+c))
##
let <- function(..., expr) {
    expr <- substitute(expr)
    dots <- list(...)
    eval(expr, dots)
}
let(a=1,b=2,expr=a+b)
## let(a=1,b=2,expr=let(c=3,expr=a+b+c))
##
let <- function(args, body) {
    f <- function() {}
    body(f) <- substitute(body)
    formals(f) <- args
    f()
}
let(list(a=1,b=2),a+b)
let(list(a=1,b=2),let(list(c=1),a+b+c)) # FAILS

let <- function(args, body) {
    bquote(function(a=1,b=3) {.(body)})
}
let(list(a=1,b=2),a+b)
let(list(a=1,b=2),let(list(c=1),a+b+c))()()


require(Rcpp)
src <- "#include <Rcpp.h>
using namespace Rcpp ;
// [[Rcpp::export]]
SEXP fn_impl( List args){
    Language call = args[0] ;
    List data = args[1] ;
    // now construct the call to eval:
    Language eval_call( \"eval\", call, data ) ;
    // and evaluate it
    return Rf_eval( eval_call, R_GlobalEnv ) ;
}"
sourceCpp(code=src)
fn_impl(list(substitute(a*2),list(a=1))) # FAILS

## stsplit/splitLexis using Rcpp
require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace std;
inline double min(double x, double y) { return (x<y ? x : y); }
inline double max(double x, double y) { return (x<y ? y : x); }
// [[Rcpp::export]]
DataFrame splitLexisCpp(NumericVector splits, NumericVector starts, NumericVector finishes) {
  vector<double> split(splits.begin(), splits.end());
  vector<int> _row, _event, _lower;
  vector<double> _t0, _t;
  int j;
  for (int i=0; i<starts.size(); ++i) {
    j = lower_bound(split.begin(), split.end(), starts[i]) - split.begin();
    if (starts[i] < split[j]) --j;
    for(; split[j] < finishes[i]; ++j) {
      _row.push_back(i+1);
      _lower.push_back(split[j]);
      _t0.push_back(max(split[j],starts[i]));
      _t.push_back(min(split[j+1],finishes[i]));
      _event.push_back(split[j]<=finishes[i] && finishes[i] < split[j+1] ? 1 : 0); // indicator for end of segment
    }
  }
  return DataFrame::create(_(".row")=wrap(_row),
                           _(".lower")=wrap(_lower),
                           _(".t0")=wrap(_t0),
                           _(".t")=wrap(_t),
                           _(".event")=wrap(_event));
}
'
sourceCpp(code=src)
splitLexis2 <- function(df,start,finish,splits) {
    start <- deparse(substitute(start))
    finish <- deparse(substitute(finish))
    stopifnot(all(df[[start]]<=df[[finish]]))
    stopifnot(max(df[[finish]])<=max(splits))
    stopifnot(min(splits)<=min(df[[start]]))
    index <- splitLexisCpp(splits,df[[start]],df[[finish]])
    cbind(df[index$.row,],index[,-1])
}
## system.time(out <- temp4 %>%
##     splitLexis2(age_start,age_finish,c(0,seq(40,90,by=5),120)) %>%
##     group_by(total_cat,.lower) %>%
##     summarise(pt = sum(.t-.t0), event=sum(event*.event)) %>%
##     mutate(rate = event/pt, age = .lower))
## NOT NEEDED: split by age group with time since last test as the primary time scale
## temp5 <- temp4 %>%
##     splitLexis2(age_start,age_finish,c(0,seq(40,90,by=5),120)) %>%
##     mutate(event = event * .event,
##            age = .lower,
##            .t0 = .t0 - age_start,
##            .t = .t - age_start)

## old versions
cureModel <- function(formula,data=NULL,formula.shape=~1,formula.cure=~1,par=NULL,...) {
    formula.scale <- formula
    rstpm2:::lhs(formula.scale) <- NULL
    y <- model.extract(model.frame(formula,data),"response")
    if (ncol(y)==2) { # right censored
        delayed <- FALSE
        time <- y[,1]
        status <- y[,2]
    }
    if (ncol(y)==3) { # left truncated and right censored
        entry <- y[,1]
        time <- y[,2]
        status <- y[,3]
        if (delayed <- any(entry>0))
            i <- which(entry>0)
    }
    X.shape <- model.matrix(formula.shape,data)
    X.scale <- model.matrix(formula.scale,data)
    X.cure <- model.matrix(formula.cure,data)
    n.shape <- ncol(X.shape)
    n.scale <- ncol(X.scale)
    n.cure <- ncol(X.cure)
    negll <-function(beta) {
        shape <- exp(as.vector(X.shape %*% beta[1:n.shape]))
        scale <- exp(as.vector(X.scale %*% beta[(1:n.scale)+n.shape]))
        cure <- rstpm2:::expit(as.vector(X.cure %*% beta[(1:n.cure)+n.scale+n.shape]))
        S <- cure+(1-cure)*pweibull(time,shape,scale,lower.tail=FALSE)
        f <- (1-cure)*dweibull(time,shape,scale)
        nll <- - sum(ifelse(status==1, log(f), log(S)))
        if (delayed) {
            S0 <- cure[i]+(1-cure[i])*pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
            nll <- nll - sum(-log(S0))
        }
        nll
    }
    if (is.null(par))
        par <- c(shape=rep(0,n.shape), # shape
                       scale=c(log(sum(status)/sum(time)),rep(0,n.scale-1)), # scale
                       cure=c(logit((1-mean(status))/2), rep(0,n.cure-1)))
    fit <- optim(par, negll, hessian=TRUE, ...)
    fit$se.fit <- sqrt(diag(solve(fit$hessian)))
    fit
}
cureModel <- function(formula,data=NULL,formula.shape=~1,formula.cure=~1,par=NULL,...) {
    formula.scale <- formula
    rstpm2:::lhs(formula.scale) <- NULL
    mf <- model.frame(formula,data)
    y <- model.extract(mf,"response")
    if (ncol(y)==2) { # right censored
        delayed <- FALSE
        time <- y[,1]
        status <- y[,2]
    }
    if (ncol(y)==3) { # left truncated and right censored
        entry <- y[,1]
        time <- y[,2]
        status <- y[,3]
        if (delayed <- any(entry>0))
            i <- which(entry>0)
    }
    formulae <- list(shape=formula.shape,
                     scale=formula.scale,
                     cure=formula.cure)
    X <- lapply(formulae,model.matrix,data)
    stopifnot(nrow(X$scale) == nrow(X$shape))
    cumn <- cumsum(sapply(X,ncol))
    index <- mapply(function(a,b) (a+1):b, c(0,cumn[-length(cumn)]),cumn,SIMPLIFY=FALSE)
    names(index) <- names(formulae)
    ll <-function(beta) {
        lp <- mapply(function(x,idx) as.vector(x %*% beta[idx]), X, index, SIMPLIFY=FALSE)
        shape <- exp(lp$shape)
        scale <- exp(lp$scale)
        cure <- rstpm2:::expit(lp$cure)
        S <- cure+(1-cure)*pweibull(time,shape,scale,lower.tail=FALSE)
        f <- (1-cure)*dweibull(time,shape,scale)
        value <- sum(ifelse(status==1, log(f), log(S)))
        if (delayed) {
            S0 <- cure[i]+(1-cure[i])*pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
            value <- value - sum(log(S0))
        }
        value
    }
    if (is.null(par))
        par <- c(shape=rep(0,n.shape), # shape
                       scale=c(log(sum(status)/sum(time)),rep(0,n.scale-1)), # scale
                       cure=c(logit((1-mean(status))/2), rep(0,n.cure-1)))
    fit <- maxLik(ll, start = par, ...)
    fit
}

require(Rcpp)
require(inline)
src <- '
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
// [[Rcpp::export]]
double logLikeCureModel(NumericVector nvtime, NumericVector nvstatus, NumericMatrix nmXshape,
NumericMatrix nmXscale, NumericMatrix nmXcure, NumericVector nvbeta) {
  double ll = 0.0;
  mat Xshape = as<mat>(wrap(nmXshape));
  mat Xscale = as<mat>(wrap(nmXscale));
  mat Xcure = as<mat>(wrap(nmXcure));
  vec time = as<vec>(wrap(nvtime));
  vec status = as<vec>(wrap(nvstatus));
  vec beta = as<vec>(wrap(nvbeta));
  int n0=Xshape.n_cols;
  int n1=n0+Xscale.n_cols;
  int n2=n1+Xcure.n_cols;
  vec shape = exp(Xshape * beta(span(0,n0-1)));
  vec scale = exp(Xscale * beta(span(n0,n1-1)));
  vec cure = 1.0/(1.0+exp(-Xcure * beta(span(n1,n2-1))));
  for (int i=0; i<nvtime.size(); ++i) {
    ll += status(i)==1.0 ?
      log(1.0-cure(i)) + R::dweibull(time(i),shape(i),scale(i),1) :
      log(cure(i)+(1.0-cure(i)) * R::pweibull(time(i),shape(i),scale(i),0,0));
  }
  return ll;
}
'
sourceCpp(code=src)
cureModel <- function(formula,data=NULL,formula.shape=~1,formula.cure=~1,par=NULL,...) {
    formula.scale <- formula
    rstpm2:::lhs(formula.scale) <- NULL
    mf <- model.frame(formula,data)
    y <- model.extract(mf,"response")
    if (ncol(y)==2) { # right censored
        delayed <- FALSE
        time <- y[,1]
        status <- y[,2]
    }
    if (ncol(y)==3) { # left truncated and right censored
        entry <- y[,1]
        time <- y[,2]
        status <- y[,3]
        if (delayed <- any(entry>0))
            i <- which(entry>0)
    }
    formulae <- list(shape=formula.shape,
                     scale=formula.scale,
                     cure=formula.cure)
    X <- lapply(formulae,model.matrix,data)
    stopifnot(nrow(X$scale) == nrow(X$shape))
    cumn <- cumsum(sapply(X,ncol))
    index <- mapply(function(a,b) (a+1):b, c(0,cumn[-length(cumn)]),cumn,SIMPLIFY=FALSE)
    names(index) <- names(formulae)
    ll <- llNew <-function(beta) {
        logLikeCureModel(time, status, X$shape,
                         X$scale, X$cure, beta)
    }
    llOld <-function(beta) {
        lp <- mapply(function(x,idx) as.vector(x %*% beta[idx]), X, index, SIMPLIFY=FALSE)
        shape <- exp(lp$shape)
        scale <- exp(lp$scale)
        cure <- rstpm2:::expit(lp$cure)
        S <- cure+(1-cure)*pweibull(time,shape,scale,lower.tail=FALSE)
        f <- (1-cure)*dweibull(time,shape,scale)
        value <- sum(ifelse(status==1, log(f), log(S)))
        if (delayed) {
            S0 <- cure[i]+(1-cure[i])*pweibull(entry[i],shape[i],scale[i],lower.tail=FALSE)
            value <- value - sum(log(S0))
        }
        value
    }
    if (is.null(par))
        par <- c(shape=rep(0,n.shape), # shape
                       scale=c(log(sum(status)/sum(time)),rep(0,n.scale-1)), # scale
                       cure=c(logit((1-mean(status))/2), rep(0,n.cure-1)))
    fit <- maxLik(ll, start = par, ...)
    fit
}


## "screening prevalence" ~= (overall screening rate)/(rescreening rate)
## = (number of screens)/pop / ((number of re-screens) / (person-time re-screened))
out <- temp4 %>%
    filter(format(start,"%Y")>="2008") %>%
    splitLexis2(age_start,age_finish,c(0,seq(40,90,by=5),120)) %>%
    group_by(.lower) %>%
    summarise(pt = sum(.t-.t0), event=sum(event*.event)) %>%
    mutate(overall_rate = event/pt, age5 = .lower)
temp2c <- temp2b %>%
    filter(year>=2008 & year<2012) %>%
    group_by(age5) %>%
    summarise(pt=sum(pt), n_tests=sum(n_tests)) %>%
    mutate(rescreen_rate=n_tests/pt)
temp2c %>%
    select(c(age5, rescreen_rate)) %>%
    inner_join(out %>% select(c(age5,overall_rate))) %>%
    mutate(prev = overall_rate / rescreen_rate)
## The rescreening rate is biased, as there is limited follow-up and declining rates for increasing duration since the previous test.

## cf. SQL
agesplits <- mutate(data.frame(start=c(0,seq(40,90,by=5))),finish=c(start[-1],120))
require(sqldf)
system.time(agg1 <- sqldf("select t1.total_cat, t2.start, t2.finish, sum(min(t1.age_finish,t2.finish) - max(t1.age_start,t2.start)) as pt, sum(event*(t1.age_finish<=t2.finish)*(t1.age_start<=t2.finish)) as event from temp4 as t1, agesplits as t2 where min(t1.age_finish,t2.finish) > max(t1.age_start,t2.start) group by t1.total_cat, t2.start, t2.finish"))

## 2D splits
require(Rcpp)
require(inline)
src <- '
#include <Rcpp.h>
#include <vector>
#include <algorithm>
using namespace Rcpp;
using namespace std;
inline double min(double x, double y) { return (x<y ? x : y); }
inline double max(double x, double y) { return (x<y ? y : x); }
// [[Rcpp::export]]
DataFrame splitLexis2DCpp(NumericVector splits0, NumericVector starts0, NumericVector finishes0,
                          NumericVector splits1, NumericVector starts1) {
  vector<double> split0(splits0.begin(), splits0.end());
  vector<double> split1(splits1.begin(), splits1.end());
  vector<int> _row, _event, _lower0, _lower1;
  vector<double> _t0, _t;
  int j, k;
  double u0, u1, offset;
  for (int i=0; i<starts0.size(); ++i) {
    offset = starts1[i] - starts0[i];
    j = lower_bound(split0.begin(), split0.end(), starts0[i]) - split0.begin();
    if (starts0[i] < split0[j]) --j;
    for(; split0[j] < finishes0[i]; ++j) {
      u0 = max(split0[j],starts0[i]) + offset;
      u1 = min(split0[j+1],finishes0[i]) + offset;
      k = lower_bound(split1.begin(), split1.end(), u0) - split1.begin();
      if (u0 < split1[k]) --k;
      for(; split1[k] < u1; ++k) {
        _row.push_back(i+1);
        _lower0.push_back(split0[j]);
        _lower1.push_back(split1[k]);
        _t0.push_back(max(split1[k],u0) - offset);
        _t.push_back(min(split1[k+1],u1) - offset);
        _event.push_back(split1[k]<=finishes0[i]+offset && finishes0[i]+offset < split1[k+1] ? 1 : 0);
      }
    }
  }
  return DataFrame::create(_(".row")=wrap(_row),
                           _(".lower1")=wrap(_lower0),
                           _(".lower2")=wrap(_lower1),
                           _(".t0")=wrap(_t0),
                           _(".t")=wrap(_t),
                           _(".event")=wrap(_event));
}
'
sourceCpp(code=src)
splitLexis2D <- function(df,start1,finish1,splits1,start2,splits2) {
    start1 <- deparse(substitute(start1))
    start2 <- deparse(substitute(start2))
    finish1 <- deparse(substitute(finish1))
    ##
    stopifnot(all(df[[start1]]<=df[[finish1]]))
    stopifnot(max(df[[finish1]])<=max(splits1))
    stopifnot(min(splits1)<=min(df[[start1]]))
    stopifnot(min(splits2)<=min(df[[start2]]))
    ##
    index <- splitLexis2DCpp(splits1,df[[start1]],df[[finish1]],splits2,df[[start2]])
    cbind(df[index$.row,],index[,-1])
}
splitLexis2D(data.frame(age0=0,age1=15,year0=1990.5),age0,age1,0:20,year0,seq(1900,2100,by=5))


require(lattice)
xyplot(rate ~ age, data=filter(out,age>=40), group=total_cat, type="l",
       auto.key=TRUE)

temp3 <- temp2 %>% filter(!is.na(total) & total>0 & age5>=50 & age5<54)
fit <- coxph(Surv(time,event)~ns(log(total),df=3),data=temp3, model=TRUE)
fit %>% summary
##
totals <- exp(with(filter(temp2,total>0),seq(min(log(total)), max(log(total)), length=301)))
X <- model.matrix(fit,data=data.frame(total=totals))
Sigma <- vcov(fit)
fitted <- as.vector(X %*% coef(fit))
se.fit <- sqrt(colSums(t(X) * (Sigma %*% t(X))))
matplot(totals,fitted+cbind(0,-1.96*se.fit,1.96*se.fit),lty=c(1,2,2),col=1,type="l",xlim=c(0,10))

require(lattice)
rmNA <- function(x) x[!is.na(x)]
d <- expand.grid(y=0:1,x=c(0,max(temp2$time/365)),group=rmNA(unique(temp2$total_cat)))
xyplot(y~x|group,
       data=d,
       panel=function(x,y,subscript) {
           plot(survfit(Surv(time/365,event)~total_cat,
                        data=temp2[subscript,],
                        subset=!is.na(total) & total>0 & age5>=70 & age5<74),
                xlab="Time from previous PSA test",
                ylab="Survival",
                col=1:4)
           })

library(survival)
library(ggplot2)
library(scales)
fortify.survfit <- function(survfit.data)
    data.frame(time = survfit.data$time,
               n.risk = survfit.data$n.risk,
               n.event = survfit.data$n.event,
               n.censor = survfit.data$n.censor,
               surv = survfit.data$surv,
               std.err = survfit.data$std.err,
               upper = survfit.data$upper,
               lower = survfit.data$lower,
               strata = rep(names(survfit.data$strata), survfit.data$strata))

d.survfit <- survfit(Surv(time/365,event)~total_cat,
                     data=temp2,
                     subset=!is.na(total) & total>0 & age5>=70 & age5<74)
ggplot(data = d.survfit) +
    geom_line(aes_string(x = 'time', y = 'surv', colour = 'strata')) +
    geom_ribbon(aes_string(x = 'time', ymin = 'lower', ymax = 'upper', fill = 'strata'), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent)
## fortified <- ggplot2::fortify(d.survfit)

do.callR <- function(args,what,...)
    do.call(what,args,...)

fortified <-
    do.call("rbind",
            lapply(seq(45,85,by=5), function(age)
                   cbind(age5=sprintf("%i-%i",age,age+4),
                         ggplot2::fortify(survfit(Surv(time/365,event)~total_cat,
                                                  data=temp2,
                                                  subset=!is.na(total) & total>0 & age5==age)))))
fortified <- seq(45,85,by=5) %>%
    lapply(function(age)
           cbind(age5=sprintf("%i-%i",age,age+4),
                 ggplot2::fortify(survfit(Surv(time/365,event)~total_cat,
                                          data=temp2,
                                          subset=!is.na(total) & total>0 & age5==age)))) %>%
    do.callR("rbind")

fortified <-
    filter(temp2, !is.na(total) & total>0 & age5>=40) %>%
    group_by(age5) %>%
    do(cbind(age5=sprintf("%i-%i",unique(.$age5),unique(.$age5)+4),
                  ggplot2::fortify(survfit(Surv(time/365,event)~total_cat,
                                           data=.)))) %>%
    mutate(strata = factor(strata, levels=c("total_cat=(0,1]", "total_cat=(1,3]","total_cat=(3,10]", "total_cat=(10,Inf]")))

ggplot(data = fortified) +
    facet_wrap( ~ age5) +
    geom_line(aes_string(x = 'time', y = 'surv', colour = 'strata')) +
    geom_ribbon(aes_string(x = 'time', ymin = 'lower', ymax = 'upper', fill = 'strata'), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent)

ggplot(data = fortified) +
    facet_wrap( ~ strata) +
    geom_line(aes_string(x = 'time', y = 'surv', colour = 'age5')) +
    geom_ribbon(aes_string(x = 'time', ymin = 'lower', ymax = 'upper', fill = 'age5'), alpha = 0.5) +
    scale_y_continuous(labels = scales::percent)





## re-screening rates - irrespective of time from last screen
## http://stackoverflow.com/questions/23026145/dplyr-how-to-number-label-data-table-by-group-number-from-group-by
group_number = (function(){i = 0; function() i <<- i+1 })()
temp4 <- temp2 %>% mutate(age_start = as.numeric((start-dob)/365), age_finish = as.numeric((finish-dob)/365)) %>% group_by(total_cat) %>% mutate(label = group_number())
temp4 <- as.data.frame(temp4)
