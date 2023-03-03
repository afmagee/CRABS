#' Stochastic exploration of congruent models.
#'
#' @param model the reference model, an object of class "CRABS"
#' @param num.samples The number of samples to be drawn
#' @param rate.type either "extinction", "speciation", "both" or "joint"
#' @param sample.speciation.rates a function that when called returns a speciation rate function
#' @param sample.extinction.rates a function that when called returns a extinction rate function
#' @param sample.joint.rates a function that when called returns a list with a speciation rate function and an extinction rate function
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
#' # Sampling extinction rates
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
#' 
#' samples
#' 
#' # Jointly sampling speciation and extinction rates
#' 
#' sample.joint.rates <- function(n) {
#'   sample.basic.models.joint(times = times, 
#'                             p.delta = model$p.delta,  
#'                             beta.param = c(0.5,0.3),  
#'                             lambda0 = l(0.0),  
#'                             mu0.median = mu(0.0))
#' }
#' 
#' joint.samples <- sample.congruence.class(model = model, 
#'                                          num.samples = 40, 
#'                                          rate.type = "joint", 
#'                                          sample.joint.rates = sample.joint.rates)
#' 
#' joint.samples
sample.congruence.class <- function(model,
                                    num.samples, 
                                    rate.type="both", 
                                    sample.speciation.rates=NULL, 
                                    sample.extinction.rates=NULL, 
                                    sample.joint.rates=NULL) {

    times <- model$times
    v_p_div <- model$p.delta(model$times)
  
    mus <- list()
    lambdas <- list()
    
    idx_lambda <- 1
    idx_mu <- 1
    
    for (i in 1:num.samples) {
      if (rate.type == "joint"){
        joint.rates = sample.joint.rates()
        lambdas[[idx_lambda]] <- joint.rates$func_lambdas
        idx_lambda <- idx_lambda + 1
        mus[[idx_mu]] <- joint.rates$func_mus
        idx_mu <- idx_mu +1
      }
      else if (rate.type == "extinction" || (rate.type == "both" && (i <= num.samples/2))) {
        mus[[idx_mu]] <- sample.extinction.rates()
        idx_mu <- idx_mu +1
      } else{
        lambdas[[idx_lambda]] <- sample.speciation.rates()
        idx_lambda <- idx_lambda + 1
      }
    }
    
    if (rate.type=="joint"){
      models <- joint.congruent.models(model, mus=mus, lambdas=lambdas)
    }else{
      models <- congruent.models(model, mus=mus, lambdas=lambdas)
    }

    return (models)
}
