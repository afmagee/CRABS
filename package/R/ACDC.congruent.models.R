ACDC.congruent.models <- function(model, mus = NULL, lambdas = NULL, keep_ref = TRUE){
  lambda0 <- model$lambda(0)
  times   <- model$times
  v_p_div <- model$p.delta( times )
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
  if (!is.null(mus)){
    
    if(length(mus) == 1 ){
      mus <- list(mus)
    }
    
    for (i in seq_along(mus)){
      cat(i,".")
      v_ext1         <- sapply(times, mus[[i]])
      
      
      v_lambda1      <- compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
      lambda_1       <- approxfun( times, v_lambda1)
      
      models1[[i]] <- congruence.class(lambda_1, mus[[i]], times, model$p.lambda, model$p.delta)
    }
    names(models1) <- paste0("model", seq_along(models1))
  }
  
  models2 <- list()
  # use lambdas to generate model
  if (!is.null(lambdas)){
    if(length(lambdas) == 1){
      lambdas <- list(lambdas)
    }
    
    
    for ( i in seq_along(lambdas)){
      lambda <- lambdas[[i]]
      if(!abs(lambda(min(times)) - model$lambda(min(times))) < 0.0001){
        stop("Initial speciation rate (at present = tips) must be equal across the congruence set.")
      }
      
      v_spec1         <- sapply(times, lambda)
      v_mu1      <- compute.extinction( v_p_div, v_spec1, delta_t )
      mu_1       <- approxfun( times, v_mu1)
      
      models2[[i]] <- congruence.class(lambda, mu_1, times, model$p.lambda, model$p.delta)
    }
    names(models2) <- paste0("model", seq_along(models2)+length(models1))
    
  }
  models <- c(models, models1, models2)
  class(models) <- c("list", "ACDCset")
  return(models)
}



#' Print method for ACDCset object
#'
#' @param models an object of class ACDCset
#' @param ... other arguments
#'
#' @return
#'
#' @examples
#' @export
print.ACDCset <- function(models, ...){
  cat("A congruent set of piecewise-linear birth-death models\n")
  cat("Knots:", length(models[[1]]$times), "\n")
  cat("Delta-tau:", models[[1]]$delta_t, "\n")
  plot.ACDCset(models)
  invisible()
}