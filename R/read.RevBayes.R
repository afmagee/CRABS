#' Title
#'
#' @param x path to log, or data frame
#' @param max_t tree height
#' @param first_n first n posterior samples
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
#' #model <- read.RevBayes("path/to/HSMRFBDP_primates.log", max_t = 73, )
#' }
read.RevBayes <- function(x,
                          max_t = 100,
                          n_samples = 20,
                          randomize = FALSE,
                          summary_type = "none",
                          extinction_prefix = "extinction_rate.",
                          speciation_prefix = "speciation_rate."){
  if(is.character(x)){
    samples <- read.table(file=x,
                          stringsAsFactors=FALSE,
                          header=TRUE)  
  }else if(is.data.frame(x)){
    samples <- x
  }

  ## Assume episodes are sorted, i.e. [1] the most recent one comes first, then [2], then [3] etc.
  speciation <- samples[, startsWith(names(samples), speciation_prefix)]
  extinction <- samples[, startsWith(names(samples), extinction_prefix)]
  
  n_epochs <- ncol(speciation)
  n_total <- nrow(speciation)
  
  times_rb <- seq(0, max_t, length.out = n_epochs)
  
  iter <- (1:n_samples) * (n_total / n_samples)
  i <- 1
  
  if (summary_type == "none"){
    pb <- txtProgressBar(min = 1, max = length(iter), style = 3)
    setTxtProgressBar(pb, 0)
    
    models <- list()

    for (it in iter){
      setTxtProgressBar(pb, i)

      lambda <- approxfun(times_rb, speciation[it,])
      mu <- approxfun(times_rb, extinction[it,])

      model <- create.model( lambda, mu, times = times_rb)
      models[[i]] <- model
      i <- i + 1
    }
    cat("\n")

    class(models) <- c("list", "ACDCposterior")
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



#' Title
#'
#' @param x 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
print.ACDCposterior <- function(x, ...){
  cat("Posterior sample\n")
  cat("Knots:", length(x[[1]]$times), "\n")
  cat("Delta-tau:", x[[1]]$delta_t, "\n")
  #p <- plot.ACDC(x, ...)
  #plot(p)
  invisible()
}
