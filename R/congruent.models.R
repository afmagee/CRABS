#' Create a set of congruent models
#'
#' @param model The reference model. An object of class "CRABS"
#' @param mus A list of extinction-rate functions
#' @param lambdas A list of speciation-rate functions
#' @param keep_ref Whether or not to keep the reference model in the congruent set
#' @param ode_solver Whether to use a numerical ODE solver to solve for lambda
#'
#' @return An object of class "CRABSset"
#' @export
#'
#' @examples
#' 
#' data(primates_ebd)
#' lambda <- approxfun(primates_ebd$time, primates_ebd$lambda)
#' mu <- approxfun(primates_ebd$time, primates_ebd$mu)
#' 
#' ## A reference model
#' times <- seq(0, max(primates_ebd$time), length.out = 500)
#' model <- create.model(lambda, mu, times = times)
#' 
#' mu1 <- lapply(c(0.5, 1.5, 3.0), function(m) function(t) m)
#' 
#' model_set1 <- congruent.models(model, mus = mu1)
#' 
#' model_set1
#' 
#' lambda0 <- lambda(0.0) ## Speciation rates must all be equal at the present
#' bs <- c(0.0, 0.01, 0.02)
#' lambda1 <- lapply(bs, function(b) function(t) lambda0 + b*t)
#' 
#' model_set2 <- congruent.models(model, lambdas = lambda1)
#' 
#' model_set2
congruent.models <- function(model, mus = NULL, lambdas = NULL, keep_ref = TRUE, ode_solver = TRUE){
  lambda0 <- model$lambda(0)
  times   <- model$times
  delta_t <- model$delta_t
  
  models <- list()
  
  if (keep_ref){
    models[[1]] <- model
    names(models) <- c("reference")
  }
  
  if(missing(lambdas)){
    if(missing(mus)){
      stop("must provide either mu(s) or lambda(s)")
    }
  }
  
  models1 <- list()
  model_idx <- 1
  # Use mus to generate model
  if (!is.null(mus) && length(mus) > 0){
    
    if(length(mus) == 1 ){
      mus <- list(mus)
    }
  
    for (i in seq_along(mus)){
      models1[[i]] <- congruent.speciation(model, mus[[i]], ode_solver = ode_solver)
      
      if (!is.null(names(mus)[[i]])){
        names(models1)[[i]] <- names(mus)[[i]]
      }else{
        names(models1)[[i]] <- paste0("model", model_idx)
        model_idx <- model_idx +1 
      }
    }
    
  }
  
  models2 <- list()
  # use lambdas to generate model
  if (!is.null(lambdas) && length(lambdas) > 0){
    if(length(lambdas) == 1){
      lambdas <- list(lambdas)
    }
    
    for ( i in seq_along(lambdas)){
      lambda <- lambdas[[i]]
      if(!abs(lambda(min(times)) - model$lambda(min(times))) < 0.0001){
        stop("Initial speciation rate (at present = tips) must be equal across the congruence set.")
      }
      
      models2[[i]] <- congruent.extinction(model, lambdas[[i]])
      if (!is.null(names(lambdas)[[i]])){
        names(models2)[[i]] <- names(lambdas)[[i]]
      }else{
        names(models2)[[i]] <- paste0("model", model_idx)
        model_idx <- model_idx +1 
      }
    }
  }
  models <- c(models, models1, models2)
  class(models) <- c("list", "CRABSset")
  return(models)
}
