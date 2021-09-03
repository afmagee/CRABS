#' Title
#'
#' @param path path to log
#' @param max_t tree height
#' @param times the time knots to be used in the piecewise-linear model
#'
#' @return a list of piecewise-linear birth-death models, each being a sample in the posterior
#' @export
#'
#' @examples
ACDC.read.RevBayes <- function(path, times = seq(from = 0, to = 73, length.out = 1000)){
  samples <- read.table(file=path, ## "output/HSMRFBDP_primates.log"
                        stringsAsFactors=FALSE,
                        header=TRUE)
  #samples <- samples[1:100,]
  
  max_t = max(times)
  
  speciation <- samples[, startsWith(names(samples), "speciation_rate.")]
  extinction <- samples[, startsWith(names(samples), "extinction_rate.")]
  
  n_epochs <- ncol(par_speciation)
  n_samples <- nrow(par_speciation)
  
  times_rb <- seq(0, max_t, length.out = n_epochs)
  pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
  setTxtProgressBar(pb, 0)
  
  models <- list()
  for (i in 1:n_samples){
    setTxtProgressBar(pb, i)
    
    lambda <- approxfun(times_rb, speciation[i,])
    mu <- approxfun(times_rb, extinction[i, ])
    
    model <- congruence.class( lambda, mu, times = times_rb)
    models[[i]] <- model
  }

  return(models)
}