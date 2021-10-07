#' Title
#'
#' @param path path to log
#' @param max_t tree height
#' @param summary_type either "none" for all the posterior samples, or "mean" or "median" for the posterior mean/median
#' @param speciation_prefix the prefix string for the speciation rate column names. Must be unique
#' @param extinction_prefix the prefix string for the extinction rate column names. Must be unique
#'
#' @return a set of ACDC models, each being a sample in the posterior
#' @export
#'
#' @examples
#' 
#' \dontrun{
#' #model <- ACDC.read.RevBayes("path/to/HSMRFBDP_primates.log", max_t = 73, )
#' }
ACDC.read.RevBayes <- function(path,
                               max_t = 100,
                               summary_type = "none",
                               extinction_prefix = "extinction_rate.",
                               speciation_prefix = "speciation_rate."){
  samples <- read.table(file=path, ## "output/HSMRFBDP_primates.log"
                        stringsAsFactors=FALSE,
                        header=TRUE)
  
  #max_t = max(times)
  
  ## Assume episodes are sorted, i.e. [1] the most recent one comes first, then [2], then [3] etc.
  speciation <- samples[, startsWith(names(samples), extinction_prefix)]
  extinction <- samples[, startsWith(names(samples), speciation_prefix)]
  
  n_epochs <- ncol(speciation)
  n_samples <- nrow(speciation)
  
  times_rb <- seq(0, max_t, length.out = n_epochs)
  
  
  if (summary_type == "none"){
    pb <- txtProgressBar(min = 0, max = n_samples, style = 3)
    setTxtProgressBar(pb, 0)
    
    models <- list()
    for (i in 1:n_samples){
      setTxtProgressBar(pb, i)

      lambda <- approxfun(times_rb, speciation[i,])
      mu <- approxfun(times_rb, extinction[i, ])

      model <- create.model( lambda, mu, times = times_rb)
      models[[i]] <- model
    }

    class(models) <- c("list", "ACDCset")
    return(models)
  }else{
    speciation_summary <- apply(speciation, 2, summary_type)
    extinction_summary <- apply(extinction, 2, summary_type)
    
    lambda <- approxfun(times_rb, speciation_summary)
    mu <- approxfun(times_rb, extinction_summary)
    
    model <- create.model( lambda, mu, times = times_rb)
    return(model)
  }
  
}