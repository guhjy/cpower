
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
#' @examples
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' one.sample(cont,means,sd=1,n=10)

#' @keywords power, contrasts, planned comparisons
#' @export

one.sample<-function(cont,means,sd=1,n=10) {

  y<-rnorm(n*length(cont),mean=means,sd=sd)
  x<-factor(rep(1:length(cont),n))
  data.frame(cbind(y,x))

}

#' Compute power based on simulation
#'
#' Simulates samples based on normal distribution and test the contrast for each sample, providing
#' power information
#'
#' @param rep number of simulations to run
#' @param cont contrast codes as a numeric vector
#' @param n per cell n
#' @param means pattern of means as numeric vector of the same length of \code{cont}
#' @param sd pooled standard deviation
#' @return a data.frame with all the simulated results of class simulations
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @examples
#'
#' ## check actual power with 10 simulations
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' pwr<-power.contr.simulate(10,cont,10,means,sd=1)
#' mean(pwr$`Pr(>|t|)`>=.05)
#'
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
class(res)<-c("simulations","data.frame")
res
}



#' Compute simulated power expected results
#'
#' Compute simulated power from power.contr.simulate() results
#'
#' @param obj a data.frame output of \code{\link[cpower]{power.contr.simulate}}
#' @seealso \code{\link[cpower]{power.contr.simulate}}
#'
#' @author Marcello Gallucci, \email{mcfanda@gmail.com}
#' @examples
#' cont<-c(-3,-1,1,3)
#' means<-c(10,12,10,12)
#' pwr<-power.contr.simulate(10,cont,10,means,sd=1)
#' summary(pwr)
#' @keywords power, contrasts, planned comparisons
#' @export

summary <- function(obj) UseMethod("summary")

summary.simulations<-function(obj) {
 rep<-dim(obj)[1]
 p<-sum((obj$`Pr(>|t|)`<.05))/rep
 t<-mean(obj$`t value`)
 c(rep=rep,power=p,ttest=t)
}
