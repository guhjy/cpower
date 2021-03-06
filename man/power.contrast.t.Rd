% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{power.contrast.t}
\alias{power.contrast.t}
\title{power function for contrasts}
\usage{
power.contrast.t(cont = NULL, d = NULL, n = NULL, power = NULL,
  scale = "g", sig.level = 0.05, type = "between",
  alternative = c("two.sided", "one.sided"))
}
\arguments{
\item{cont}{contrast codes as a numeric vector}

\item{d}{effect size d}

\item{n}{number of observations per cell}

\item{power}{power of test (1 minus Type II error probability)}

\item{scale}{used to scale the d index. It can be "g", "z" or a numeric value}

\item{sig.level}{significance level (Type I error probability)}

\item{type}{string specifying the type of t test. Can be abbreviated.}

\item{alternative}{one- or two-sided test. Can be abbreviated}
}
\value{
Object of class '"power.htest"', a list of the arguments (including the computed one) augmented with 'method' and 'note' elements.
}
\description{
computes power parameters for t-test associated with a custom contrast
}
\details{
Exactly one of the parameters 'd','n','power'
     and 'sig.level' must be passed as NULL, and that parameter is determined from the others.

     The parameter \code{scale} controls the method used to scale the effect size d.
     \enumerate{
    \item    \code{scale="g"} assumes scaling by dividing 2*d by the sum of absolute coefficients
    \item    \code{scale="z"} assumes scaling by dividing d by the square-root of the sum of squares of the contrast weigths
    \item    \code{numeric} any constant that multiplies the unscaled d to obtain the scaled d
   }
}
\author{
Marcello Gallucci, \email{mcfanda@gmail.com}
}
\keyword{comparisons}
\keyword{contrasts,}
\keyword{planned}
\keyword{power,}
