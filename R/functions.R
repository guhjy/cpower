#' build a custom contrasts and k-2 orthogonal ones
#'
#' Creates a set of contrast codes that are orthogonal to a given one with the aim of
#' testing the given contrast in a linear model
#'
#' @param cont Character string indicating a set of colors.
#'
#' @return a data.frame with k-1 contrast codes, with k=length(cont).
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons

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


test.contr<-function(data,yname,xname,cont,debug=FALSE) {
  con<-contr.custom(cont)
  data[,xname]<-factor(data[,xname])
  contrasts(data[,xname])<-con
  form<-as.formula(paste(yname,"~",xname))
  ss<-summary(lm(form,data=data))
  if (debug)
    print(ss)
  round(ss$coefficients[2,],5)
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
#'     \item    \code{scale="q"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons



power.contrast.t<-function(cont=NULL,d=NULL,n=NULL,power=NULL,scale="g",sig.level=0.05) {

  if (scale=="g") {
    w<-sum(abs(cont))/2
    d0=d*w
  } else if (scale=="q") {
    d0=d*sqrt(sum(cont^2))
  } else {
    d0=d/scale
  }
  k=length(cont)
  ss<-sum(cont^2)
  alpha<-sig.level
  p.body <- quote({
    ncp<-sqrt(n/ss)*d0
    np<-(n-1)*k
    crit<-qt(alpha/2,np,lower.tail = F)
    pt(crit, np, ncp = ncp, lower.tail = FALSE) +  pt(-crit, np, ncp = ncp, lower.tail = TRUE)
  })
  if (is.null(power)) {
    power <- eval(p.body)
    return(power)
  }
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
    return(n)
  }
}

#' power function for contrasts
#'
#' computes power parameters for F-test associated with a custom contrast
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
#'     \item    \code{scale="q"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#'
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons


power.contrast.f<-function(cont=NULL,d=NULL,n=NULL,power=NULL,scale="g",sig.level=0.05) {

  if (scale=="g") {
    w<-sum(abs(cont))/2
    d0=d*w
  } else if (scale=="q") {
    d0=d*sqrt(sum(cont^2))
  } else {
    d0=d/scale
  }
  k=length(cont)
  ss<-sum(cont^2)
  alpha<-sig.level
  p.body <- quote({
    ncp<-((n/ss)*d0^2)
    np<-(n-1)*k
    crit<-qf(sig.level, 1, np, lower.tail = FALSE)
    pf(crit,1, np, ncp,lower.tail = FALSE)
  })
  if (is.null(power)) {
    power <- eval(p.body)
    return(power)
  }
  if (is.null(n)) {
    n <- uniroot(function(n) eval(p.body) - power, c(2, 1e+05))$root
    return(n)
  }
}



