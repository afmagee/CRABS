#' Title
#'
#' @param posterior a list of ACDC model objects
#' @inheritParams sample.congruence.class
#'
#' @inherit sample.congruence.class return
#' @export
#'
#' @examples 
#' data(primates_ebd_log)
#' 
#' posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 20)
#' 
#' extinction_rate_samples <- function(){
#' res <- sample.basic.models(
#'     num.epochs = 100,
#'     rate0.median = 0.1,
#'     model = "MRF",
#'     max.rate = 1.0)
#'   return(res)
#' }
#' 
#' samples <- sample.congruence.class.posterior(posterior, 
#'                                              num.samples = 20,
#'                                              rate.type = "extinction",
#'                                              sample.extinction.rates = extinction_rate_samples)
#' 
#' print(samples)
sample.congruence.class.posterior <- function(posterior,
                                              num.samples, 
                                              rate.type="extinction", 
                                              sample.extinction.rates=NULL){
  if(rate.type == "speciation" || rate.type == "both"){
    stop("speciation rate samples not currently implemented")
  }
  
  pb <- txtProgressBar(min = 1, max = length(posterior), style = 3)  
  res <- list()
  
  for(i in seq_along(posterior)){
    cg <- sample.congruence.class(posterior[[i]],
                                  num.samples = num.samples,
                                  rate.type = "extinction",
                                  sample.extinction.rates = sample.extinction.rates)
    res[[i]] <- cg
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  names(res) <- paste0("posterior", seq_along(res))
  class(res) <- c("ACDCsets", "list")
  return(res)
}

#' Title
#'
#' @param posterior a list of ACDC objects, each one representing a sample from the posterior
#' @inheritParams summarize.trends
#' 
#' @return a ggplot object
#' @export
#'
#' @examples
#' data(primates_ebd_log)
#' 
#' posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 20)
#' 
#' extinction_rate_samples <- function(){
#' res <- sample.basic.models(
#'     num.epochs = 100,
#'     rate0.median = 0.1,
#'     model = "MRF",
#'     max.rate = 1.0)
#'   return(res)
#' }
#' 
#' samples <- sample.congruence.class.posterior(posterior, 
#'                                              num.samples = 20,
#'                                              rate.type = "extinction",
#'                                              sample.extinction.rates = extinction_rate_samples)
#' 
#' 
#' summarize.posterior(samples, threshold = 0.05)
summarize.posterior <- function(posterior,
                                threshold = 0.01,
                                rate_name = "lambda",
                                return_data = FALSE,
                                rm_singleton = FALSE,
                                relative_deltas = FALSE){
  times <- posterior[[1]][[1]]$times
  max_t <- max(times)
  n_epochs <- length(times)
  
  res <- list()
  for (i in seq_along(posterior)){
    df <- summarize.trends(posterior[[i]], 
                           threshold = threshold, 
                           return_data = TRUE)[["heatmap_data"]]
    df["posterior"] <- paste0("sample",i)
    res[[i]] <- df
  }
  plotdata <- bind_rows(res)
  
  if(return_data){
    return(df)
  }
  
  k <- sum(sapply(posterior, length))
  #scaleyformat <- function(x) sprintf("%.1f", x/k)

  p1 <- ggplot(plotdata, aes(x = time, fill = direction)) +
    geom_histogram(aes(y = stat(count / sum(count))), binwidth = max_t/n_epochs) +
    theme_classic() +
    scale_fill_manual(values = c("purple","white", "#7fbf7b"), labels = direction_labels) +
    theme(legend.position = c(0.3, 0.3),
          legend.title = element_blank()) +
    scale_x_reverse() +
    labs(y = "model coverage", x = "time before present")
  return(p1)
}



#' print.ACDCsets
#'
#' @param x a list of (congruent) ACDC sets
#' @param ... additional parameters
#'
#' @return nothing
#' @export
#'
#' @examples
#' data(primates_ebd_log)
#' 
#' posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 20)
#' 
#' extinction_rate_samples <- function(){
#' res <- sample.basic.models(
#'     num.epochs = 100,
#'     rate0.median = 0.1,
#'     model = "MRF",
#'     max.rate = 1.0)
#'   return(res)
#' }
#' 
#' samples <- sample.congruence.class.posterior(posterior, 
#'                                              num.samples = 20,
#'                                              rate.type = "extinction",
#'                                              sample.extinction.rates = extinction_rate_samples)
#'                                              
#' print(samples)
print.ACDCsets <- function(x, ...){
  cat("A group of ", length(x), "ACDC sets.\n")
  cat("Knots:", length(x[[1]][[1]]$times), "\n")
  cat("Delta-tau:", x[[1]][[1]]$delta_t, "\n")
  invisible()
}