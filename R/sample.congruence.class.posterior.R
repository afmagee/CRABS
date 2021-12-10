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