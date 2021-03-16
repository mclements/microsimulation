#' summary method for a SummaryReport object
#'
#' @param object SummaryReport object
#' @param ... other arguments
#' @return a list of class summary.SummaryReport with components:
#' \describe{
#' \item{n}{Number of simulations}
#' \item{indivip}{boolean with whether individual values were retained}
#' \item{utilityDiscountRate}{discount rate for utilities/QALYs}
#' \item{costDiscountRate}{discount rate for costs}
#' \item{QALE}{Quality-adjusted life expectancy (discounted)}
#' \item{LE}{Life expectancy (not discounted)}
#' \item{ECosts}{Life-time expected costs (discounted)}
#' \item{se.QALE}{standard error for QALE}
#' \item{se.Ecosts}{standard error Ecosts}
#' }
#' @rdname SummaryReport
#' @export
summary.SummaryReport = function(object, ...)
    with(object,
         structure(list(n = n,
                        indivp = indivp,
                        utilityDiscountRate = utilityDiscountRate,
                        costDiscountRate = costDiscountRate,
                        QALE = sum(ut$utility)/n,
                        LE = sum(pt$pt)/n,
                        Ecosts = sum(costs$cost)/n,
                        se.QALE = sd(indiv$utilities)/sqrt(n),
                        se.Ecosts = sd(indiv$costs)/sqrt(n)),
                   class="summary.SummaryReport"))

#' Print summary from SummaryReport object
#'
#' @param x summary.SummaryReport object
#' @param ... other arguments passed to print
#' @rdname SummaryReport
#' @export
print.summary.SummaryReport <- function(x,...)
    with(x,
         print(c("n"=n,"Utility discount rate"=utilityDiscountRate,"Cost discount rate"=costDiscountRate,"Cost"=Ecosts,"(se)"=se.Ecosts,"QALYs"=QALE,
                 "(se)"=se.QALE),
               ...))

#' Print SummaryReport object
#'
#' At present, this passes the object to summary and then prints
#' 
#' @param x SummaryReport object
#' @param ... other arguments passed to print
#' @rdname SummaryReport
#' @export
print.SummaryReport <- function(x,...)
    print(summary(x),...)

#' Row bind a set of SummaryReport objects
#'
#' @param ... a set of SummaryReport objects
#' @return a SummaryReport object
#' @rdname SummaryReport
#' @export
rbind.SummaryReport <- function(...) {
    objects = list(...)
    stopifnot(all(sapply(objects, function(obj) obj$param$utilityDiscountRate)==
                  objects[[1]]$param$utilityDiscountRate))
    stopifnot(all(sapply(objects, function(obj) obj$param$costDiscountRate)==
                  objects[[1]]$param$costDiscountRate))
    newobject = objects[[1]]
    newobject$n = sum(sapply(objects, "[[", "n"))
    for (name in c("pt","ut","events","prev","costs","indiv"))
        newobject[[name]] = do.call(rbind,lapply(objects, "[[", name))
    newobject
}

#' ascii output from a SummaryReport
#'
#' @param object a SummaryReport object
#' @param include.rownames logical for whether to include rownames (default=FALSE)
#' @param include.colnames logical for whether to include colnames (default=TRUE)
#' @param header logical for whether to include the header (default=TRUE)
#' @param digits vector of the number of digits to use for each column
#' @param ... other arguments to pass to ascii
#' @return ascii object
#' @rdname SummaryReport
#' @export
ascii.SummaryReport <- function(object,include.rownames=FALSE,include.colnames=TRUE,header=TRUE,
                                digits=c(0,3,2,2,4,4),...) {
    if (requireNamespace("ascii")) {
        with(summary(object),
             ascii(c("n"=n,"Discount rate"=discountRate,"Cost"=Ecosts,"(se)"=se.Ecosts,"QALYs"=QALE,
                     "(se)"=se.QALE),
                   include.rownames, include.colnames, header=header, digits=digits, ...))
    } else stop("ascii package not available")
}

#' ICER for two SummaryReport objects
#'
#' @param object1 SummaryReport object (reference)
#' @param object2 SummaryReport object
#' @param ... other arguments (not currently used)
#' @return a list of type ICER.SummaryReport with components:
#' \describe{
#' \item{n}{number of simulations}
#' \item{utilityDiscountRate}{Discount rate for the utilities/QALE}
#' \item{costDiscountRate}{Discount rate for the costs}
#' \item{s1}{summary for object1}
#' \item{s2}{summary for object2}
#' \item{dQALE}{QALE for object2 minus QALE for object1}
#' \item{dCosts}{Costs for object2 minus costs for object1}
#' \item{ICER}{change of costs divided by change in QALEs}
#' \item{se.dQALE}{standard error for dQALE}
#' \item{se.dCosts}{standard error for dCosts}
#' }
#' @rdname SummaryReport
#' @export
ICER.SummaryReport = function(object1, object2, ...) {
    stopifnot(object1$n == object2$n)
    stopifnot(object1$utilityDiscountRate == object2$utilityDiscountRate)
    stopifnot(object1$costDiscountRate == object2$costDiscountRate)
    stopifnot(object1$indivp == object2$indivp)
    s1 = summary(object1)
    s2 = summary(object2)
    dQALE = s2$QALE - s1$QALE
    dCosts = s2$Ecosts - s1$Ecosts
    structure(list(n=object1$n,
                   utilityDiscountRate=object1$utilityDiscountRate,
                   costDiscountRate=object1$costDiscountRate,
                   s1=s1, s2=s2,
                   dQALE=dQALE, dCosts=dCosts, ICER=dCosts/dQALE,
                   se.dQALE = sd(object2$indiv$utilities-object1$indiv$utilities)/sqrt(object1$n),
                   se.dCosts = sd(object2$indiv$costs-object1$indiv$costs)/sqrt(object1$n)),
              class="ICER.SummaryReport")
}

#' ascii output from a ICER.SummaryReport object
#'
#' @param object an ICER.SummaryReport object
#' @param include.rownames logical for whether to include rownames (default=FALSE)
#' @param include.colnames logical for whether to include colnames (default=TRUE)
#' @param header logical for whether to include the header (default=TRUE)
#' @param digits vector of the number of digits to use for each column
#' @param rownames rownames for output
#' @param colnames colnames for output
#' @param tgroup tgroup arg passed to ascii
#' @param n.tgroup arg passed to ascii
#' @param ... other arguments to pass to ascii
#' @return ascii object
#' @rdname SummaryReport
#' @export
ascii.ICER.SummaryReport <-
    function(object,include.rownames=TRUE,include.colnames=TRUE,header=TRUE,
             digits=c(1,1,3,3,1,1,3,3,1),
             rownames=c("Reference","Treatment"),
             colnames=c("Costs","(se)","QALYs","(se)","Costs","(se)","QALYs","(se)","ICER"),
             tgroup=c("Total","Incremental"),n.tgroup=c(4,5),...) {
        if (requireNamespace("ascii")) {
            m <- with(object,
                      matrix(c(s1$Ecosts,s1$se.Ecosts,s1$QALE,s1$se.QALE,NA,NA,NA,NA,NA,
                               s2$Ecosts,s2$se.Ecosts,s2$QALE,s2$se.QALE,dCosts,se.dCosts,
                               dQALE, se.dQALE, ICER),2,byrow=TRUE))
            dimnames(m) = list(rownames,colnames)
            ascii(m,include.rownames,include.colnames,header=header,digits=digits,
                  tgroup=tgroup, n.tgroup=n.tgroup, ...)
            } else stop("ascii package not available")
    }

