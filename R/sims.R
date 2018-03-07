
#' Return a sample with the specified pattern of means and standard deviation
#'
#' Simulates samples based on normal distribution and test the contrast for each sample, providing
#' power information
#'
#' @param cont contrast codes as a numeric vector
#' @param means pattern of means as numeric vector of the same length of \code{cont}
#' @param sd pooled standard deviation
#' @param n per cell n
#' @return a data.frame with all the simulated results of class simulations
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export


one.sample<-function(cont,means,sd=1,n=100) {

  y<-rnorm(n*length(cont),mean=means,sd=sd)
  x<-factor(rep(1:length(cont),n))
  data<-data.frame(cbind(y,x))


}

#' Compute power based on simulations
#'
#' Simulates samples based on normal distribution and test the contrast for each sample, providing
#' power information
#'
#' @param cont contrast codes as a numeric vector
#' @param means pattern of means as numeric vector of the same length of \code{cont}
#' @param sd pooled standard deviation
#' @param n per cell n
#' @return a data.frame with all the simulated results of class simulations
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export

power.contr.simulate<-function(rep,cont,n,means,sd=1) {

res<-replicate(rep, {
  y<-rnorm(n*length(cont),mean=means,sd=sd)
  x<-factor(rep(1:length(cont),n))
  data<-data.frame(cbind(y,x))
  test.contr(data,"y","x",cont)
},simplify = T)
res<-as.data.frame(t(rbind(res)))
class(res)<-c("simulation","data.frame")
res
}

#' Compute simulated power expected results
#'
#' Compute simulated power from power.contr.simulate() results
#'
#' @param obj a data.frame outout of \code{\link[cpower]{power.contr.simulate}}
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @seealso \code{\link[cpower]{cpower}}
#' @keywords power, contrasts, planned comparisons
#' @export

summary.simulations<-function(obj) {
 rep<-dim(obj)[1]
 p<-sum((obj$`Pr(>|t|)`<.05))/rep
 t<-mean(obj$`t value`)
 c(rep=rep,power=p,ttest=t)
}
