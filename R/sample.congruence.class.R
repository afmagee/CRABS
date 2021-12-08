#' Stochastic exploration of congruent models.
#'
#' @param model the reference model, an object of class "ACDC"
#' @param num.samples The pulled diversification rate function (measured in time before present).
#' @param rate.type either "extinction", "speciation", or "both"
#' @param sample.speciation.rates a function that when called returns a vector representing the speciation rate trajectory
#' @param sample.extinction.rates a function that when called returns a vector representing the extinction rate trajectory
#' @return A named list with congruent rates.
#' @export
#' @examples
#' #TODO

sample.congruence.class <- function(model,
                                    num.samples, 
                                    rate.type="both", 
                                    sample.speciation.rates=NULL, 
                                    sample.extinction.rates=NULL) {
  
    times <- model$times
    v_p_div <- model$p.delta(times)
  
    mus <- list()
    lambdas <- list()
    
    pb <- txtProgressBar(min = 0, max = num.samples, style = 3)
    setTxtProgressBar(pb, 0)

    idx_lambda <- 1
    idx_mu <- 1
    
    for (i in 1:num.samples) {
      if (rate.type == "extinction" || (rate.type == "both" && (i <= num.samples/2))) {
        this_mu     <- sample.extinction.rates()
        func_ext1   <- approxfun(times, this_mu)
        
        mus[[idx_mu]] <- func_ext1
        idx_mu <- idx_mu +1
      } else {
        this_lambda <- sample.speciation.rates()
        func_spec1  <- approxfun(times,this_lambda)
        
        lambdas[[idx_lambda]] <- func_spec1
        idx_lambda <- idx_lambda + 1
      }
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    models <- congruent.models(model, mus = mus, lambdas = lambdas)
    
    #cat("Sampled ",n_rates_drawn," rate trajectories to get ",num.samples," valid extinction rates (",round(100*num.samples/n_rates_drawn,2),"%)\n",sep="")
    cat("Sampled ", num.samples, " rate trajectories.\n")
    
    return (models)
}
