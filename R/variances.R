#' Compute (partial) eta-squared for contrast
#'
#' Compute the eta-squared for a contrast based on a model already estimated with the effect of
#' "xname" af a factor.
#'
#' @param cont  vector of contrast weights.
#' @param xname Character string indicating the name of the factor defining groups.
#' @param model a \code{\link[stats]{lm}} object with \code{xname} as independent variable.
#'
#' @return A numeric value for eta-squared.
#' @description
#' Return the proportion of variance explained by the contrast over the variance not explained by other effects.
#' The result is the partial eta-squared unless, there's only one possible contrast in the model, the is the eta-squared
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link[cpower]{eta2.contr.d}}
#' @examples
#' ### sim some data
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' y<-rep(means,1000)+rnorm(4000,0,2)
#' x<-factor(rep(1:4,1000))
#' ### compute eta-squared
#' mod<-lm(y~x)
#' eta2.contr(cont,xname = "x",model = mod)
#'
#' #### check from the d coefficient
#' d<-d.contr(cont,x=x,y=y)
#' eta2.contr.d(cont = cont,d=d,scale = "g")
#' @keywords power, contrasts, planned comparisons, eta-squared
#' @export

eta2.contr<-function(cont,xname,model) {
  .data<-model$model
  .data[,xname]<-factor(.data[,xname])
  n<-dim(.data)[1]
  contrasts(.data[,xname])<-contr.custom(cont)
  new<-update(model,data=.data)
  ss<-summary(new)
  co<-ss$coefficients
  err<-model$df.residual*(ss$sigma)^2
  .name<-paste(xname,"cont",sep="")
  param<-co[.name,][1]*sum(cont^2)
  ssc<-n*param^2/((sum(cont^2))*length(cont))
  ssc/(ssc+err)
}

#' Compute (partial) eta-squared for contrast d
#'
#'
#' @param cont  vector of contrast weights.
#' @param d d coefficient for given contrast.
#' @param scale method to scale d.
#' @return eta-squared.
#' @details
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="q"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#' @description
#' Return the proportion of variance explained by the contrast over the variance not explained by other effects.
#' The result is the partial eta-squared unless, there's only one possible contrast in the model, the is the eta-squared
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link[cpower]{eta2.contr}}
#' @examples
#' ### sim some data
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' y<-rep(means,1000)+rnorm(4000,0,2)
#' x<-factor(rep(1:4,1000))
#' ### compute eta-squared
#' mod<-lm(y~x)
#' eta2.contr(cont,xname = "x",model = mod)
#'
#' #### check from the d coefficient
#' d<-d.contr(cont,x=x,y=y)
#' eta2.contr.d(cont = cont,d=d,scale = "g")

#' @keywords power, contrasts, planned comparisons
#' @export
eta2.contr.d<-function(cont,d,scale="g") {
  d0<-.tod0(cont,d,scale)
  dq<-d0/sqrt(sum(cont^2))
  dq^2/(dq^2+length(cont))
}


#' Compute contrast d from (partial) eta-squared
#'
#'
#' @param cont  vector of contrast weights.
#' @param eta2 eta-squared coefficient for given contrast.
#' @param scale method to scale d.
#' @return eta-squared.
#'  @details
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="q"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @export
d.contr.eta2<-function(cont,eta2,scale="g") {
  dz<-sqrt((eta2*length(cont))/(1-eta2))
  if (scale=="g") {
    w<-sum(abs(cont))/2
    d=dz*sqrt(sum(cont^2))/w
  } else if (scale=="z") {
    d=dz
  } else {
    d=dz*sqrt(sum(cont^2))/scale
  }
d
}


#' Compute contrast f from d
#'
#'
#' @param cont  vector of contrast weights.
#' @param d d coefficient for given contrast.
#' @param scale method to scale d.
#' @return f.
#' @details
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @export
f.contr.d<-function(cont,d,scale="g") {
  d0<-.tod0(cont,d,scale)
  dz<-.fromd0(cont,d0,"z")
  abs(dz)/sqrt(length(cont))
}


#' Compute contrast d from f
#'
#'
#' @param cont  vector of contrast weights.
#' @param f f coefficient for given contrast.
#' @param scale method to scale d.
#' @return f.
#' @details
#'      The parameter \code{scale} controls the method used to scale the effect size d.
#'      \enumerate{
#'     \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
#'     \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the coefficients
#'     \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
#'    }
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @keywords power, contrasts, planned comparisons
#' @export
d.contr.f<-function(cont,f,scale="g") {
  dz<-f*sqrt(length(cont))
  d0<-.tod0(cont,dz,"z")
  .fromd0(cont,d0,scale)
}

