#' Compute confidence interval for contrast d
#'
#' Compute the confidence interval based on the Steiger and Fouladi (1997) method
#'
#' @param cont  vector of contrast weights.
#' @param d standardized mean differences.
#' @param scale method to scale d.
#' @param n cell size (for each group)
#' @param conf.level width of the confidence interval, default=.95
#' @details
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'      The confidence intervals are computed based on selected function in \code{MBESS}:
#'      in particular \code{\link[MBESS]{conf.limits.nct}} and \code{\link[MBESS]{ci.sc}}. Results are scaled to the
#'      standardized effect size required.

#' @return confidence intervals of class "conf.intervals".
#' @description Return the confidence intervals computed by scaling the confidence intervals of the noncentrality parameter.
#' The confidence intervals are computed based on selected function in \code{MBESS}.
#' @references
#' Kelley, K. (2007). Confidence intervals for standardized effect sizes: Theory, application, and implementation. Journal of Statistical Software, 20 (8), 1–24.
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @examples
#'
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' d<-d.contr(cont,means = means,sd=2,scale = "g")
#' ci.contr(cont,d=d,scale = "g",n=100)
#'
#' @export
#'
#'


ci.contr<-function(cont,d,n,scale="g",conf.level=.95) {
  sub<-ifelse(is.numeric(scale),"c",scale)

  k<-length(cont)
  d0<-.tod0(cont,d,scale=scale)
  ncp<-d0*(sqrt(n/sum(cont^2)))
  ci<-ci.ncpt(ncp = ncp,k = k,n = n,conf.level = conf.level )
  ci0<-lapply(ci,function(x) x*sqrt(sum(cont^2)/n))
  sci<-lapply(ci0,.fromd0,cont=cont,scale=scale)
  structure(list(lower=as.numeric(sci[1]),
                 param=as.numeric(sci[2]),
                 upper=as.numeric(sci[3])), class="conf.intervals",conf.level=conf.level,df=n*(k-1),param=paste("d",sub,sep=""))

}



#' Compute confidence interval for noncentral t
#'
#' Compute the confidence interval based on the Steiger and Fouladi (1997) method for the noncentrality parameter of t
#'
#' @param ncp  noncentrality parameter
#' @param k number of groups.
#' @param n cell size (for each group)
#' @param conf.level width of the confidence interval, default=.95

#' @return confidence intervals.
#' @description Return the confidence intervals computed by scaling the confidence intervals of the noncentrality parameter.
#' @references
#'   Steiger, J. H., & Fouladi, R. T. (1997). Noncentrality interval estimation and the evaluation of statistical models. In L. L. E. Harlow, S. A. Mulaik, & J. H. Steiger (Eds.), What if there were no significance tests?: Classic edition. Lawrence Erlbaum Associates Publishers.
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @export
#'
#'


ci.ncpt<-function(ncp,k,n,conf.level=.95) {

  .find.ncp <- function(val)
  {
    qt(p=alpha, df=k*(n-1), ncp=val, lower.tail = lowertail) - ncp
  }
  alpha=(1-conf.level)/2
  suppressWarnings({
    lowertail<-FALSE
    low<-try(uniroot(f=.find.ncp,interval=c(-37, 37)),silent = T)
    if (class(low)=="try-error") {
      low=-Inf
    }
    lowertail<-TRUE
    up<-try(uniroot(f=.find.ncp,interval=c(-37, 37)),silent = T)
    if (class(up)=="try-error") {
      low=Inf
    }
  })
  structure(list(lower=as.numeric(low[[1]]),
            param=ncp,
            upper=as.numeric(up[[1]])), class="conf.intervals",conf.level=conf.level,df=n*(k-1),param="Ncp")
}


#' Nice printing for confidence intervals
#'
#'
#' @param obj  an object of class "conf.intervals"

#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link[cpower]{ci.contr}}
#' @keywords power, contrasts, planned comparisons
#' @export
#'
#'
print <- function(obj) UseMethod("print")


#' Nice printing for confidence intervals
#'
#'
#' @param obj  an object of class "conf.intervals"

#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link[cpower]{ci.contr}}
#' @keywords power, contrasts, planned comparisons
#' @export print.conf.intervals
#'
#'

print.conf.intervals<-function(obj) {
  alpha<-attr(obj,"conf.level")
  center<-attr(obj,"param")
  df<-attr(obj,"df")
  cat(round(alpha*100),"% Confidence interval (df=",df,")\n\n")
  res<-as.table(c(obj$lower,obj$param,obj$upper))
  names(res)<-c("Lower",center,"Upper")
  print(res)
  cat("\n")
}
