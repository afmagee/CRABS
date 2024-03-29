% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.congruence.class.posterior.R
\name{sample.congruence.class.posterior}
\alias{sample.congruence.class.posterior}
\title{Stochastic exploration of congruent models for all samples in the posterior}
\usage{
sample.congruence.class.posterior(
  posterior,
  num.samples,
  rate.type = "extinction",
  mu0.equal = FALSE,
  rate0 = NULL,
  ...
)
}
\arguments{
\item{posterior}{a list of CRABS model objects}

\item{num.samples}{The number of samples to be drawn}

\item{rate.type}{either "extinction", "speciation", "both" or "joint"}

\item{mu0.equal}{whether to propose alternative mu starting at mu0 equal to the posterior sample. default to FALSE}

\item{rate0}{rate0 allows the user to fix the extinction rate at the present to a single value. defaults to NULL, for drawing it randomly}

\item{...}{
  Arguments passed on to \code{\link[=sample.basic.models]{sample.basic.models}}
  \describe{
    \item{\code{times}}{the time knots}
    \item{\code{model}}{"MRF" for pure MRF model, otherwise MRF has a trend of type "exponential","linear", or "episodic<n>"}
    \item{\code{direction}}{"increase" or "decrease" (measured in past to present)}
    \item{\code{noisy}}{If FALSE, no MRF noise is added to the trajectory}
    \item{\code{MRF.type}}{"HSMRF" or "GMRF", type for stochastic noise.}
    \item{\code{monotonic}}{Whether the curve should be forced to always move in one direction.}
    \item{\code{fc.mean}}{Determines the average amount of change when drawing from the model.}
    \item{\code{rate0.median}}{When not specified, rate at present is drawn from a lognormal distribution with this median.}
    \item{\code{rate0.logsd}}{When not specified, rate at present is drawn from a lognormal distribution with this sd}
    \item{\code{mrf.sd.scale}}{scale the sd of the mrf process up or down. defaults to 1.0}
    \item{\code{min.rate}}{The minimum rate (rescaling fone after after drawing rates).}
    \item{\code{max.rate}}{The maximum rate (rescaling fone after after drawing rates).}
  }}
}
\value{
A named list with congruent rates.
}
\description{
This function takes a posterior sample as input: a list of CRABS objects. 
It will then iterate over the samples, and for each posterior sample it will
sample from the posterior class. It will sample using the \code{\link{sample.basic.models}}
function, and all additional parameters are passed to \code{\link{sample.basic.models}}.
}
\examples{
data(primates_ebd_log)

posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 10)

samples <- sample.congruence.class.posterior(posterior, 
                                             num.samples = 5,
                                             rate.type = "extinction",
                                             rate0.median = 0.1,
                                             model = "MRF",
                                             max.rate = 1.0)

print(samples)
}
