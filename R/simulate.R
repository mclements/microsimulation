#' Simulate event times from a survreg object
#' @param object survreg object
#' @param nsim number of simulations per row in newdata
#' @param seed random number seed
#' @param newdata data-frame for defining the covariates for the simulations. Required.
#' @param t0 delayed entry time. Defaults to NULL (which assumes that t0=0)
#' @param ... other arguments (not currently used)
#' @return vector of event times with nsim repeats per row in newdata
#' @importFrom stats simulate predict runif
#' @importFrom survival survreg.distributions
#' @rdname simulate
#' @export
#' @examples
#' library(survival)
#' fit <- survreg(Surv(time, status) ~ ph.ecog + age + sex + strata(sex),
#'                data = lung)
#' nd = transform(expand.grid(ph.ecog=0:1, sex=1:2), age=60)
#' simulate(fit, seed=1002, newdata=nd)
#' simulate(fit, seed=1002, newdata=nd, t0=500)
simulate.survreg = function(object, nsim=1, seed=NULL, newdata, t0=NULL, ...) {
    stopifnot(inherits(object, "survreg"),
              is.list(newdata))
    if (!is.null(seed)) set.seed(seed)
    lp = predict(object, newdata=newdata, type="lp")
    lp = lp[rep(1:length(lp), each=nsim)]
    n = length(lp)
    if (!is.null(strata <- attr(object$terms, "specials")$strata)) {
        scale = object$scale[eval(attr(object$terms,"variables")[[strata+1]], newdata)]
        scale = rep(scale, each=nsim)
    }
    else scale = rep(object$scale,n)
    if (is.character(object$dist)) 
        dd <- survival::survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
        trans <- function(x) x
        itrans <- function(x) x
    }
    else {
        trans <- dd$trans
        itrans <- dd$itrans
    }
    if (!is.null(dd$dist)) 
        dd <- survival::survreg.distributions[[dd$dist]]
    if (!is.null(t0)) {
        stopifnot(length(t0) %in% c(1, length(newdata[[1]])))
        if (length(t0)==1) t0=rep(t0,n)
        else t0 = rep(t0,each=nsim)
        F0 = dd$density((trans(t0)-lp)/scale,object$parm)[,1]
    } else F0 = rep(0,n)
    itrans(lp+scale*dd$quantile(1-(runif(n, F0, 1)-F0), object$parm))
}
