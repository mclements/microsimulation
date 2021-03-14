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

print.summary.SummaryReport <- function(x,...)
    with(x,
         print(c("n"=n,"Utility discount rate"=utilityDiscountRate,"Cost discount rate"=costDiscountRate,"Cost"=Ecosts,"(se)"=se.Ecosts,"QALYs"=QALE,
                 "(se)"=se.QALE),
               ...))

print.SummaryReport <- function(x,...)
    print(summary(x),...)

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

ascii.SummaryReport <- function(object,include.rownames=FALSE,include.colnames=TRUE,header=TRUE,
                                digits=c(0,3,2,2,4,4),...) {
    if (requireNamespace("ascii")) {
        with(summary(object),
             ascii(c("n"=n,"Discount rate"=discountRate,"Cost"=Ecosts,"(se)"=se.Ecosts,"QALYs"=QALE,
                     "(se)"=se.QALE),
                   include.rownames, include.colnames, header=header, digits=digits, ...))
    } else stop("ascii package not available")
}

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

