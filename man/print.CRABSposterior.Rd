% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.RevBayes.R
\name{print.CRABSposterior}
\alias{print.CRABSposterior}
\title{Title}
\usage{
\method{print}{CRABSposterior}(x, ...)
}
\arguments{
\item{x}{a list of CRABS objects}

\item{...}{additional parameters}
}
\value{
nothing
}
\description{
Title
}
\examples{
data(primates_ebd_log)
posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 20)
print(posterior)
}
