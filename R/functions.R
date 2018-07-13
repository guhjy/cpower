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
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the contrast weigths
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @examples
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' d.contr(cont,means = means,sd=2,scale = "g")
#' ### different scaling
#' d.contr(cont,means = means,sd=2,scale = "z")
#' ### raw data
#' y<-rep(means,1000)+rnorm(4000,0,2)
#' x<-factor(rep(1:4,1000))
#' d.contr(cont,y = y,x = x)
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
    sd<-sqrt(sum(v*n)/df)
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


#' Build a custom contrast and k-2 orthogonal ones
#'
#' Creates a set of contrast codes that are orthogonal to a given one with the aim of
#' testing the given contrast in a linear model
#'
#' @param cont vector of contrast weight (numeric).
#'
#' @return a data.frame with k-1 contrast codes, with k=length(cont). It can be used to assing contrast weights using \code{\link{contrasts}}
#' @examples
#' ### sim some data
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' y<-rep(means,1000)+rnorm(4000,0,2)
#' x<-factor(rep(1:4,1000))
#' ### assign contrast weights
#' contrasts(x)<-contr.custom(cont)
#' contrasts(x)
#' cor(contr.custom(cont))
#' summary(lm(y~x))

#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
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
#' Perform a t-test for custom contrast
#' @param cont contrast codes as a numeric vector with length=k
#' @param d effect size index to test
#' @param n total sample size required if d is provided
#' @param scale scaling method used to compute d. Required if d is provided
#' @param data data.frame containing the variables to be analysed
#' @param yname name (string) of the dependent variable
#' @param xname name (string) of the independent variable with k-groups, with k being length(cont)
#'
#' @return t-test for the contrast. The \code{Estimate} is the contrast value estimate from the data and the \code{Std. Error} is expressed in the same
#'         scale of the estimate. It is assumed that the cells have the same size.
#' @details The test can be performed either on raw data on a d coefficient. When \code{d} is provided, also the total sample size \code{n}
#'          and the scaling method \code{scale} should be provided.
#'
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the contrast weigths
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }

#'
#' @examples
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' y<-rep(means,1000)+rnorm(4000,0,10)
#' x<-rep(1:4,1000)
#' dat<-as.data.frame(cbind(y,x))
#' dat$x<-factor(dat$x)
#'
#' test.contr(data = dat,cont=cont,yname = "y",xname = "x")
#'
#' # check the contrast value
#' observed<-tapply(dat$y,dat$x,mean)
#' observed%*%cont
#'
#' # check the t-test and p-value
#' contrasts(dat$x)<-contr.custom(cont)
#' summary(lm(y~x,data=dat))

#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons, t-test
#' @export


test.contr<-function(cont,d=NULL,n=NULL,scale=NULL,y=NULL,x=NULL) {
  if (is.null(d)) {
     if (is.null(y) | is.null(y))
        stop("Either d or y and x vectors should be provided")
      con<-contr.custom(cont)
      x<-factor(x)
      contrasts(x)<-con
      form<-as.formula(y~x)
      model<-lm(form)
      ss<-summary(model)
      res<-ss$coefficients[2,]
      res[1]<-res[1]*sum(cont^2)
      res[2]<-res[2]*sum(cont^2)
     return(res)
  }
  if (is.null(y) & is.null(x)) {
    if (is.null(d))
      stop("Either d or y and x vectors should be provided")
    if (is.null(n))
      stop("When d is provided, also the total sample size n is required")
    if (is.null(scale))
      stop("When d is provided, also the scaling method is required")
   k<-length(cont)
   d0<-.tod0(cont,d,scale)
   dz<-.fromd0(cont,d0,"z")
   ttest<-sqrt(n/k)*dz
   p<-2*pt(-abs(ttest),df=n-k)
   return(c(d=d,"t value"=ttest,"Pr(>|t|)" =p))
  }

 stop("Either d or x and y vectors should be provided")

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
#' @param  scale used to scale the d index. It can be "g", "z" or a numeric value
#' @param type string specifying the type of t test. Can be abbreviated.
#' @param alternative one- or two-sided test. Can be abbreviated
#' @return Object of class '"power.htest"', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
#' @details Exactly one of the parameters 'd','n','power'
#'      and 'sig.level' must be passed as NULL, and that parameter is determined from the others.
#'
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the contrast weigths
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @export

power.contrast.t<-function(cont=NULL,d=NULL,n=NULL,power=NULL,scale="g",sig.level=0.05,type="between",alternative=c("two.sided","one.sided")) {

  tol = .Machine$double.eps^0.25
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
    } else if (scale=="z") {
      d0=d*sqrt(sum(cont^2))
    } else {
      d0=d/scale
    }
  }

  k=length(cont)
  ss<-sum(cont^2)

  p.body <-
        quote({
          .ncp<-sqrt(n/ss)*d0
          np<-(n-1)*k
          crit<-qt(sig.level/tside,np,lower.tail = F)
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

  if (is.null(sig.level))
       sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10), tol = tol, extendInt = "yes")$root

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
    w<-2/sum(abs(cont))
    return(d/w)
  }
  if (scale=="z")
      return(d*sqrt(sum(cont^2)))
  return(d0/scale)
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
