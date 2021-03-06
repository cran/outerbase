% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{loglik_std}
\alias{loglik_std}
\alias{Rcpp_loglik_std-class}
\alias{Rcpp_loglik_std}
\title{Gaussian errors}
\arguments{
\item{om}{an \code{\link{outermod}} instance to be referred to}

\item{terms}{a matrix of \code{terms}, must have as many columns as dims in 
\code{om}}

\item{y}{a vector of observations}

\item{x}{a matrix of predictors, must have as many columns as dims in 
\code{om} and the same number of rows as \code{y}}
}
\value{
no returns, this is a class which contains methods
}
\description{
\preformatted{
loglik = new(loglik_std, om, terms, y, x)
}
This is a standard model which has the form
\deqn{y = \langle \phi(x), \theta \rangle + \varepsilon, \varepsilon \sim 
N(0,\sigma^2)}
where \eqn{\phi(x)} is the basis, \eqn{\theta} is the coefficient vector,
\eqn{\varepsilon} is an unseen noise vector.
The parameter vector is of length 1 where 
\code{para} \eqn{= \log(\sigma)}.  It is a slower (sometimes) 
version of \code{\link{loglik_gauss}} but allows for complete marginal 
inference.
}
\seealso{
base class: \code{\link{lpdf}}
}
