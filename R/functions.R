#' Cohen's d for contrasts
#'
#' Compute Cohen's d for contrasts with different scaling functions
#'
#' @param cont contrast codes as a numeric vector
#' @param means pattern of means as numeric vector of the same length of \code{cont}
#' @param sd pooled standard deviation
#' @param y  dependent variable
#' @param x  independent variable with k-groups, with k being length(cont)
#' @param scale scaling method of the d index
#' @return t-test for the contrast.
#' @details data to compute the d index can be either actual data or a vector of means and a standard deviation
#'          actual data are specified with the \code{y} and \code{x} parameters, means and standard deviation
#'          with the \code{means} , and \code{sd} parameters. The parameter \code{means} and \code{sd} are prevalent, thus
#'          if \code{means},\code{y} and \code{x} are provided, but not \code{sd}, y and x are used to computed the pooled standard deviation.
#'          If \code{means} and \code{sd} are provided, \code{y} and \code{x} are ignored.
#'
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export


d.contr<-function(cont,means=NULL,sd=NULL,y=NULL,x=NULL,scale="g") {
  if (is.null(cont)) stop("contrast weights must be provided")
  .badinput<-"data should be provided either as a vector of means and a pooled sd or as x and y variables"
  if (all(is.null(means),is.null(y),is.null(x)))
       stop(.badinput)
  if (is.null(means)) {
    if (is.null(y) | is.null(x))
       stop(.badinput)
    .x<-factor(x)
    means<-tapply(y,.x,mean, na.rm=T)
  }
  if (is.null(sd)) {
    if (is.null(y) | is.null(x))
      stop(.badinput)
    .x<-factor(x)
    n<-tapply(y,.x,length)
    n<-n-1
    v<-tapply(y,.x,var)
    df=sum(n)
    sd<-sum(v*n)/df
  }

  d<-.method0(cont,means,sd)
  scaled=FALSE
  if (scale=="g") {
     d<-.methodg(cont,d)
     scaled=TRUE
  }
  if (scale=="z") {
    d<-.methodq(cont,d)
    scaled=TRUE
  }
  if (is.numeric(scale)) {
     d<-d*scale
     scaled=TRUE
  }
  if (!scaled)
    warning("The index was not scaled")
   return(as.numeric(d))

}


#' build a custom contrast and k-2 orthogonal ones
#'
#' Creates a set of contrast codes that are orthogonal to a given one with the aim of
#' testing the given contrast in a linear model
#'
#' @param cont vector of contrast weight (numeric).
#'
#' @return a data.frame with k-1 contrast codes, with k=length(cont).
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export


contr.custom<-function(cont) {
  k<-length(cont)
  nc<-k-1
  if (nc==1)
    return(cont)
  data<-data.frame(cont)
  base<-matrix(runif(k*(nc-1),0,10000),ncol=nc-1)
  base<-data.frame(cbind(cont,base))
  for (i in 1:(nc-1)) {
    nnames<-names(base)
    form1<-paste(nnames[1:i],collapse = "+")
    form<-as.formula(paste(nnames[i+1],form1,sep = "~"))
    data[,i+1]<-residuals(lm(form,data=base))
  }
  as.matrix(data)
}

#' t-test for a custom contrast
#'
#' t-test for custom contrast
#'
#' @param data data.frame containing the variables to be analysed
#' @param yname name (string) of the dependent variable
#' @param xname name (string) of the independent variable with k-groups, with k being length(cont)
#' @param cont contrast codes as a numeric vector with length=k
#' @param debug if TRUE, prints out the full model results. If FALSE (default) prints only the contrast results
#'
#' @return t-test for the contrast.
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export


test.contr<-function(data,yname,xname,cont,debug=FALSE) {
  con<-contr.custom(cont)
  data[,xname]<-factor(data[,xname])
  contrasts(data[,xname])<-con
  form<-as.formula(paste(yname,"~",xname))
  model<-lm(form,data=data)
  ss<-summary(model)
  ss$coefficients[2,]*sum(cont^2)
}

print.test.contr<-function(obj) {
print(summary.lm(obj))
}


#' power function for contrasts
#'
#' computes power parameters for t-test associated with a custom contrast
#'
#' @param cont contrast codes as a numeric vector
#' @param d effect size d
#' @param n number of observations per cell
#' @param sig.level significance level (Type I error probability)
#' @param power power of test (1 minus Type II error probability)
#' @param  method used to scale the d index
#' @return Object of class '"power.htest"', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
#' @details Exactly one of the parameters 'd','n','power'
#'      and 'sig.level' must be passed as NULL, and that parameter is determined from the others.
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export

power.contrast.t<-function(cont=NULL,d=NULL,n=NULL,power=NULL,scale="g",sig.level=0.05,type="between",alternative=c("two.sided","one.sided")) {

  .ncp<-NULL
  NOTE<-""
  if (!is.numeric(cont))
     stop("'cont' must be a numeric vector")
  if (sum(sapply(list(n, d, sd, power, sig.level), is.null)) != 1)
    stop("exactly one of 'n', 'delta', 'sd', 'power', and 'sig.level' must be NULL")
  if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > sig.level | sig.level > 1))
          stop("'sig.level' must be numeric in [0, 1]")
  type <- match.arg(type)
  alternative <- match.arg(alternative)
  #tsample <- switch(type, between = k, repeated = 1)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)

  if (tside == 2 && !is.null(d))
    d <- abs(d)

  if (!is.null(d)) {
    if (scale=="g") {
      w<-sum(abs(cont))/2
      d0=d*w
    } else if (scale=="q") {
      d0=d*sqrt(sum(cont^2))
    } else {
      d0=d/scale
    }
  }

  k=length(cont)
  ss<-sum(cont^2)
  alpha<-sig.level


  p.body <-
        quote({
          .ncp<-sqrt(n/ss)*d0
          np<-(n-1)*k
          crit<-qt(alpha/tside,np,lower.tail = F)
          pt(crit, np, ncp = .ncp, lower.tail = FALSE) +  pt(-crit, np, ncp = .ncp, lower.tail = TRUE)
      })

  if (is.null(power)) {
    power <- eval(p.body)
  }
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
  }
  if (is.null(d))
    d <- uniroot(function(d0) eval(p.body) - power,
                      c(1e-07, 1e+07), extendInt = "upX")$root

  NOTE <- switch(type, paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs",
                 two.sample = "n is number in *each* group", NULL)
  structure(list(n = n, d = d, k=k, sig.level = sig.level,
                 power = power, alternative = alternative, note = NOTE,
                 scale = scale,method="Contrast t-test power calculation"), class = "power.htest")
}


.method0<-function(cont,means,sd) {
  (cont%*%means)/sd
}

.methodg<-function(cont,d0) {
  g<-2/sum(abs(cont))
  g*d0
}

.methodq<-function(cont,d0) {
  q<-1/sqrt(sum(cont^2))
  q*d0
}

.tod0<-function(cont,d,scale) {
  if (scale=="g") {
    w<-sum(abs(cont))/2
    d0=d*w
  } else if (scale=="z") {
    d0=d*sqrt(sum(cont^2))
  } else {
    d0=d/scale
  }
d0
}
.fromd0<-function(cont,d0,scale) {
  if (scale=="g") {
    w<-2/sum(abs(cont))
    d=d0*w
  } else if (scale=="z") {
    d=d0/sqrt(sum(cont^2))
  } else {
    d=d0*scale
  }
  d
}
