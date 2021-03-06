% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{gethyp}
\alias{gethyp}
\title{Get the hyperparameters}
\arguments{
\item{om}{an \code{\link{outermod}} instance}
}
\value{
a vector of parameters
}
\description{
\preformatted{
hyp = gethyp(om)
}
Gets the current hyperparameters from an \code{\link{outermod}} instance. It 
formats them in a way that makes reading in \code{R} easier.
}
\examples{
om = new(outermod)
setcovfs(om, c("mat25", "mat25", "mat25"))
hyp = gethyp(om)
print(hyp)
}
\seealso{
\code{\link{outermod}}
}
