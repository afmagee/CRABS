#' Sample custom functions through time.
#'
#' @param num.epochs The number of rates to draw
#' @param lambda0 The rate at present
#' @param rsample Function to sample next rate
#' @param rsample0 Function to sample rate at present
#' @param autocorrelated Should rates be autocorrelated?
#' @return Sampled rate vector
#' @export
#' @example 
#' #TODO
sample.rates <- function(num.epochs, lambda0=NULL, rsample=NULL, rsample0=NULL, autocorrelated=FALSE) {
  
  N_SAMPLES = ifelse( is.null(lambda0), num.epochs+1, num.epochs )
  if ( autocorrelated == FALSE ) {
    # we draw a bunch of iid samples
    new_rates = rsample(N_SAMPLES)
    if( is.null(lambda0) == FALSE ) {
      new_rates = c(lambda0, new_rates)
    }
  } else {
    # we draw autocorrelated rates
    if ( is.null(lambda0) == FALSE ) {
      new_rates = c( lambda0 )
    } else if ( is.null(rsample0) == FALSE ) {
      new_rates = c( rsample0() )
    }
    for ( i in 1:num.epochs ) {
      new_rates[i+1] = rsample(new_rates[i])
    }
  }
  
  return (new_rates)
}


#' Samples simple increase/decrease models through time with noise.
#'
#' @param num.epochs The number of rates to draw
#' @param rate0 The rate at present, otherwise drawn randomly.
#' @param model "MRF" for pure MRF model, otherwise MRF has a trend of type "exponential","linear", or "episodic<n>"
#' @param MRF.type "HSMRF" or "GMRF", type for stochastic noise.
#' @param direction "increase" or "decrease" (measured in past to present)
#' @param monotonic Whether the curve should be forced to always move in one direction.
#' @param noisy If FALSE, no MRF noise is added to the trajectory
#' @param fc.mean Determines the average amount of change when drawing from the model.
#' @param rate0.median When not specified, rate at present is drawn from a lognormal distribution with this median.
#' @param rate0.sd When not specified, rate at present is drawn from a lognormal distribution with this sd
#' @param minimum.rate The minimum rate (rescaling fone after after drawing rates).
#' @param maximum.rate The maximum rate (rescaling fone after after drawing rates).
#' @return Extinction rate at all changepoints
#' @export
#' @example 
#' #TODO
sample.basic.models <- function(num.epochs, rate0=NULL, model="exponential", MRF.type="HSMRF", direction="decrease", noisy=TRUE, fc.mean=3, rate0.median=0.1, rate0.sd=1.17481, monotonic=FALSE, min.rate=0, max.rate=10) {
  # recover()
  
  s <- 1
  if ( direction == "increase" ) {
    s <- -1
  } else if ( direction != "decrease" ) {
    stop("Invalid \"direction\"")
  }
  
  # rate at present
  x0 <- min.rate - 10
  while ( x0 < min.rate || x0 > max.rate ) {
    if ( is.null(rate0) ) {
      x0 <- rlnorm(1,log(rate0.median),rate0.sd)
    } else {
      x0 <- rate0
    }
  }
  
  fc_mean_adj <- fc.mean - 1
  fc_rate <- 1.25
  fc_shape <- fc_mean_adj * fc_rate
  fc <- 1 + rgamma(1,fc_shape,fc_rate)
  # cat(fc,"\n")
  
  # Deterministic component
  delta_deterministic <- rep(0,num.epochs-1)
  x <- numeric(num.epochs)
  x[1] <- x0
  if ( model == "exponential" ) {
    delta_deterministic <- rep(s*log(fc)/(num.epochs-1),num.epochs-1)
    x[2:num.epochs] <- x[1] * exp(cumsum(delta_deterministic))
  } else if ( model == "linear" ) {
    if (s == 1) {
      delta_deterministic <- rep(((x0*fc)-x0)/(num.epochs-1),num.epochs-1)
    } else {
      delta_deterministic <- rep(((x0/fc)-x0)/(num.epochs-1),num.epochs-1)
    }
    x[2:num.epochs] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("episodic",model) ) {
    njumps <- as.numeric(gsub("episodic","",model)) - 1
    if (njumps < 1) {
      stop("Too few episodes in episodic model")
    }
    if ( njumps == 1 ) {
      delta_deterministic[sample.int(num.epochs-1,1)] <- fc
    } else {
      delta_deterministic[sample.int(num.epochs-1,njumps)] <- fc * rdirichlet(njumps,1)
    }
    x[2:num.epochs] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("MRF",model) ) {
    x <- rep(x0,num.epochs)
  } else {
    stop("Invalid \"model\"")
  }
  
  # Add noise
  if ( noisy ) {
    # Get random component of rate changes
    zeta <- 0
    delta_stochastic <- rep(Inf,num.epochs-1)
    noise <- rep(Inf,num.epochs)
    # Avoid infinities
    while ( any(!is.finite(noise)) ) {
      if ( MRF.type == "HSMRF"  || MRF.type == "HSRF") {
        zeta <- get.hsmrf.global.scale(num.epochs)
        gamma <- abs(rcauchy(1,0,1))
        gamma <- min(gamma,1000) # avoid numerical instability
        sigma <- abs(rcauchy(num.epochs-1,0,1))
        delta_stochastic <- rnorm(num.epochs-1,0,sigma*gamma*zeta)
      } else if ( MRF.type == "GMRF" ) {
        zeta <- get.gmrf.global.scale(num.epochs)
        gamma <- min(gamma,1000) # avoid numerical instability
        gamma <- abs(rcauchy(1,0,1))
        delta_stochastic <- rnorm(num.epochs-1,0,gamma*zeta)
      } else {
        stop("Invalid \"MRF.type\"")
      }
      
      if ( monotonic ) {
        delta_stochastic <- s * abs(delta_stochastic)
      }
      
      noise <- c(1,exp(cumsum(delta_stochastic)))
    }
    x <- x * noise
  }
  
  # Make sure minimum and maximum are respected
  if ( any(x > max.rate) ) {
    adj <- (max(x) - x0)/(max.rate - x0)
    x[x > x0] <- x0 + (x[x > x0] - x0) / adj
  }
  if ( any(x < min.rate) ) {
    adj <- abs(min(x) - x0)/abs(min.rate - x0)
    x[x < x0] <- x0 + (x[x < x0] - x0) / adj
  }
  
  return (x)
}
