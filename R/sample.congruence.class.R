#' Stochastic exploration of congruent models.
#'
#' @param model the reference model, an object of class "CRABS"
#' @param num.samples The pulled diversification rate function (measured in time before present).
#' @param rate.type either "extinction", "speciation", or "both"
#' @param sample.speciation.rates a function that when called returns a speciation rate function
#' @param sample.extinction.rates a function that when called returns a extinction rate function
#' @return A named list with congruent rates.
#' @export
#' @examples
#' data("primates_ebd")
#' 
#' l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' times <- primates_ebd[["time"]]
#' 
#' model <- create.model(l, mu, primates_ebd[["time"]])
#' 
#' extinction_rate_samples <- function(){
#'    res <- sample.basic.models(times = times, 
#'                               rate0 = 0.05, 
#'                               model = "MRF", 
#'                               MRF.type = "HSMRF", 
#'                               fc.mean = 2.0, 
#'                               min.rate = 0.0, 
#'                               max.rate = 1.0)
#'    return(res)
#' } 
#' 
#' samples <- sample.congruence.class(model, 
#'                                    num.samples = 8,
#'                                    rate.type = "extinction",
#'                                    sample.extinction.rates = extinction_rate_samples)
sample.congruence.class <- function(model,
                                    num.samples, 
                                    rate.type="both", 
                                    sample.speciation.rates=NULL, 
                                    sample.extinction.rates=NULL) {

    times <- model$times
    v_p_div <- model$p.delta(model$times)
  
    mus <- list()
    lambdas <- list()
    
    idx_lambda <- 1
    idx_mu <- 1
    
    for (i in 1:num.samples) {
      if (rate.type == "extinction" || (rate.type == "both" && (i <= num.samples/2))) {
        mus[[idx_mu]] <- sample.extinction.rates()
        idx_mu <- idx_mu +1
      } else {
        lambdas[[idx_lambda]] <- sample.speciation.rates()
        idx_lambda <- idx_lambda + 1
      }
      
    }
    models <- congruent.models(model, mus = mus, lambdas = lambdas)

    return (models)
}
