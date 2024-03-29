% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample.rates.R
\name{sample.rates}
\alias{sample.rates}
\title{Sample custom functions through time.}
\usage{
sample.rates(
  times,
  lambda0 = NULL,
  rsample = NULL,
  rsample0 = NULL,
  autocorrelated = FALSE
)
}
\arguments{
\item{times}{the time knots}

\item{lambda0}{The rate at present}

\item{rsample}{Function to sample next rate}

\item{rsample0}{Function to sample rate at present}

\item{autocorrelated}{Should rates be autocorrelated?}
}
\value{
Sampled rate vector
}
\description{
Sample custom functions through time.
}
\examples{
data("primates_ebd")

l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
times <- primates_ebd[["time"]]

model <- create.model(l, mu, times)

rsample <- function(n) runif(n, min = 0.0, max = 0.9)
mu <- sample.rates(times, 0.5, rsample = rsample)


model_set <- congruent.models(model, mus = mu)

model_set
}
