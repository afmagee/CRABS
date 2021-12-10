#' Create a set of congruent models
#'
#' @param model The reference model. An object of class "ACDC"
#' @param mus A list of extinction-rate functions
#' @param lambdas A list of speciation-rate functions
#' @param keep_ref Whether or not to keep the reference model in the congruent set
#'
#' @return An object of class "ACDCset"
#' @export
#'
#' @examples
#' 
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' ## A reference model
#' model <- create.model(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' mu1 <- lapply(c(0.5, 1.5, 3.0), function(m) function(t) m)
#' 
#' model_set1 <- congruent.models(model, mus = mu1)
#' 
#' model_set1
#' 
#' lambda0 = lambda(0.0) ## Speciation rates must all be equal at the present
#' bs <- c(-0.1, 0.0, 0.1)
#' lambda1 <- lapply(bs, function(b) function(t) lambda0 + b*t)
#' 
#' model_set2 <- congruent.models(model, lambdas = lambda1)
#' 
#' model_set2
congruent.models <- function(model, mus = NULL, lambdas = NULL, keep_ref = TRUE){
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
  # Use mus to generate model
  if (!is.null(mus) && length(mus) > 0){
    
    if(length(mus) == 1 ){
      mus <- list(mus)
    }
  
    for (i in seq_along(mus)){
      models1[[i]] <- congruent.speciation(model, mus[[i]])
    }
    names(models1) <- paste0("model", seq_along(models1))
    
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
    }
    names(models2) <- paste0("model", seq_along(models2)+length(models1))
    
  }
  models <- c(models, models1, models2)
  class(models) <- c("list", "ACDCset")
  return(models)
}