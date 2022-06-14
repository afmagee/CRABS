#' read RevBayes log file
#'
#' @param x path to log, or data frame
#' @param n_times number of time knots
#' @param max_t tree height
#' @param n_samples first n posterior samples
#' @param summary_type either "none" for all the posterior samples, or "mean" or "median" for the posterior mean/median
#' @param speciation_prefix the prefix string for the speciation rate column names. Must be unique
#' @param extinction_prefix the prefix string for the extinction rate column names. Must be unique
#' @usage read.RevBayes(x, n_times, max_t = 100, n_samples = 20, summary_type = "none", 
#' extinction_prefix = "extinction_rate.", speciation_prefix = "speciation_rate.")
#'
#' @return a set of CRABS models, each being a sample in the posterior
#' @export 
#'
#' @examples
#' data(primates_ebd_log)
#' posterior <- read.RevBayes(primates_ebd_log, n_times = 500, max_t = 65, n_samples = 20)
read.RevBayes <- function(x,
                          n_times,
                          max_t = 100,
                          n_samples = 20,
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
  if(missing("n_times")){
    n_times <- length(times_rb)
  }
  
  iter <- floor(seq(1, nrow(samples), length.out = n_samples))
  i <- 1
  
  if (summary_type == "none"){
    pb <- txtProgressBar(min = 1, max = length(iter), style = 3)
    setTxtProgressBar(pb, 0)
    
    models <- list()

    for (it in iter){
      setTxtProgressBar(pb, i)
      
      if (any(is.na(speciation[it,])) || any(is.na(extinction[it,])) || any(speciation[it,] < 0) || any(extinction[it,] < 0)){
        warning(paste("Posterior sample", it," containing negative or NA rate values. skipping."))
      }else{
        lambda <- approxfun(times_rb, speciation[it,])
        mu <- approxfun(times_rb, extinction[it,])
        times <- seq(0, max(times_rb), length.out = n_times)
        
        model <- create.model( lambda, mu, times = times)
        models[[i]] <- model
        i <- i + 1
      }

      
    }
    close(pb)
    cat("\n")

    names(models) <- paste0("posterior", seq_along(models))
    class(models) <- c("list", "CRABSposterior")
    return(models)
  }else{
    speciation_summary <- apply(speciation, 2, summary_type)
    extinction_summary <- apply(extinction, 2, summary_type)
    
    lambda <- approxfun(times_rb, speciation_summary)
    mu <- approxfun(times_rb, extinction_summary)
    times <- seq(0, max(times_rb), length.out = n_times)
    
    model <- create.model( lambda, mu, times = times)
    return(model)
  }
  
}



#' Title
#'
#' @param x a list of CRABS objects
#' @param ... additional parameters
#'
#' @return nothing
#' @export
#'
#' @examples
#' data(primates_ebd_log)
#' posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 20)
#' print(posterior)
print.CRABSposterior <- function(x, ...){
  cat("Posterior sample\n")
  cat("Knots:", length(x[[1]]$times), "\n")
  cat("Delta-tau:", x[[1]]$delta_t, "\n")
  #p <- plot.CRABS(x, ...)
  #plot(p)
  invisible()
}
