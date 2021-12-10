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