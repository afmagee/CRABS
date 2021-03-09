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
#' @param direction "increase" or "decrease" (measured in past to present)
#' @param noisy If FALSE, no MRF noise is added to the trajectory
#' @param MRF.type "HSMRF" or "GMRF", type for stochastic noise.
#' @param monotonic Whether the curve should be forced to always move in one direction.
#' @param fc.mean Determines the average amount of change when drawing from the model.
#' @param rate0.median When not specified, rate at present is drawn from a lognormal distribution with this median.
#' @param rate0.logsd When not specified, rate at present is drawn from a lognormal distribution with this sd
#' @param minimum.rate The minimum rate (rescaling fone after after drawing rates).
#' @param maximum.rate The maximum rate (rescaling fone after after drawing rates).
#' @return Speciation or extinction rate at a number of timepoints.
#' @export
#' @example 
#' #TODO
sample.basic.models <- function(num.epochs=100, rate0=NULL, model="exponential", direction="decrease", noisy=TRUE, MRF.type="HSMRF", monotonic=FALSE, fc.mean=3, rate0.median=0.1, rate0.logsd=1.17481, min.rate=0, max.rate=10) {
  # recover()
  s <- 1
  if ( direction == "increase" ) {
    s <- -1
  } else if ( direction != "decrease" ) {
    stop("Invalid \"direction\"")
  }
  
  # We use rejection sampling to find a model that fits within minimum and maximum rates
  # To speed up sampling, we break this into three rejection sampling steps, handling separately the:
  #    1) rate at present
  #    2) fold change (given rate at present)
  #    3) MRF noise
  
  # rate at present
  x0 <- min.rate - 10
  if ( !is.null(rate0) ) {
    if ( x0 < min.rate || x0 > max.rate ) {
      stop("User-defined rate0 is outside [min.rate,max.rate].")
    }
  } else {
    while ( x0 < min.rate || x0 > max.rate ) {
      x0 <- rlnorm(1,log(rate0.median),rate0.logsd)
    }
  }
  
  # draw a fold change
  fc_mean_adj <- fc.mean - 1
  fc_rate <- 1.25
  fc_shape <- fc_mean_adj * fc_rate
  fc <- Inf
  while ( x0 * fc < min.rate || x0 * fc > max.rate ) {
    fc <- 1 + rgamma(1,fc_shape,fc_rate)
  }
  # cat(fc,"\n")
  # cat(fc * x0,"\n")
  
  # Deterministic component of trajectory
  delta_deterministic <- rep(0,num.epochs)
  x <- numeric(num.epochs+1)
  x[1] <- x0
  if ( model == "exponential" ) {
    delta_deterministic <- rep(s*log(fc)/(num.epochs),num.epochs)
    x[2:(num.epochs+1)] <- x[1] * exp(cumsum(delta_deterministic))
  } else if ( model == "linear" ) {
    if (s == 1) {
      delta_deterministic <- rep(((x0*fc)-x0)/(num.epochs),num.epochs)
    } else {
      delta_deterministic <- rep(((x0/fc)-x0)/(num.epochs),num.epochs)
    }
    x[2:(num.epochs+1)] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("episodic",model) ) {
    njumps <- as.numeric(gsub("episodic","",model)) - 1
    if (njumps < 1) {
      stop("Too few episodes in episodic model")
    }
    if ( njumps == 1 ) {
      delta_deterministic[sample.int(num.epochs,1)] <- (fc * x0) - x0
    } else {
      delta_deterministic[sample.int(num.epochs,njumps)] <- (fc * x0) - x0 * rdirichlet(njumps,1)
    }
    x[2:(num.epochs+1)] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("MRF",model) ) {
    x <- rep(x0,num.epochs+1)
  } else {
    stop("Invalid \"model\"")
  }
  
  rates <- x
  
  # Add noise
  if ( noisy ) {
    found_valid_model <- FALSE
    while ( !found_valid_model ) {
      # Get random component of rate changes
      zeta <- 0
      delta_stochastic <- rep(Inf,num.epochs)
      noise <- rep(Inf,num.epochs)
      # Avoid infinities, this loop should rarely trigger
      while ( any(!is.finite(noise)) ) {
        if ( MRF.type == "HSMRF"  || MRF.type == "HSRF") {
          zeta <- get.hsmrf.global.scale(num.epochs+1)
          gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
          sigma <- abs(rcauchy(num.epochs,0,1))
          delta_stochastic <- rnorm(num.epochs,0,sigma*gamma*zeta)
        } else if ( MRF.type == "GMRF" ) {
          zeta <- get.gmrf.global.scale(num.epochs+1)
          gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
          delta_stochastic <- rnorm(num.epochs,0,gamma*zeta)
        } else {
          stop("Invalid \"MRF.type\"")
        }
        
        if ( monotonic ) {
          delta_stochastic <- s * abs(delta_stochastic)
        }
        
        noise <- c(1,exp(cumsum(delta_stochastic)))
      }
      
      x_prop <- x * noise
      
      if ( all(x_prop <= max.rate) && all(x_prop >= min.rate) ) {
        rates <- x * noise
        found_valid_model <- TRUE
      }
    }
  }

  return (rates)
}

#' Samples curves like the given curve by bootstrapping.
#'
#' @param num.epochs The discretization fineness
#' @param func_rate Function for getting rates at time t.
#' @param max.t The maximum time (before present) to consider rates. 
#' @param keep.rate0 If TRUE, rate at present is fixed, otherwise chosen as func_rate at a random time.
#' @param reverse.time If TRUE, resampled curve is "backwards" from real curve, such that increases become decreases and vice-versa. Incompatible with keep.rate0=TRUE
#' @param invert If TRUE, returns a value proportional to 1/rates, allowing increases to become decreases while preserving rate at present.
#' @param replace If FALSE, only permutations of the curve are considered, if TRUE curves will be more variable.
#' @param block.size If replace=TRUE, larger values allow resampling to preserve larger-scale features of the rate trajectory.
#' @return Rate at all changepoints
#' @export
#' @example 
#' #TODO
bootstrap.rates <- function(num.epochs, func_rate, max.t, keep.rate0=TRUE, reverse.time=FALSE, invert=FALSE, replace=TRUE, block.size=1) {
  # recover()
  
  rates <- func_rate(seq(0,max.t,length.out=num.epochs+1))
  deltas <- rates[-1] - rates[-(num.epochs+1)]
  
  x0 <- ifelse(keep.rate0,rates[1],rates[sample.int(num.epochs+1,1)])
  
  resampled_deltas <- rep(0,num.epochs)
  
  if ( block.size != 1 ) {
    if ( block.size < 1 || block.size > num.epochs) {
      stop("Invalid block.size, need 1 <= block.size <= num.epochs")
    }
    n_blocks <- ceiling(num.epochs/block.size)
    n_full_blocks <- floor(num.epochs/block.size)
    blocks <- list()
    if ( replace ) {
      block_sizes <- c()
      if ( n_blocks * block.size == num.epochs ) {
        block_sizes <- rep(block.size,n_blocks)
      } else {
        leftover <- num.epochs - ((n_full_blocks)*block.size) + 1
        block_sizes <- c(rep(block.size,n_full_blocks),leftover)
      }
      blocks <- lapply(1:n_blocks,function(i){
        start <- sample.int(num.epochs-block_sizes[i]+1,1)
        return((start):(start+block_sizes[i]-1))
      })
    } else {
      starts <- (0:(n_full_blocks-1))*block.size
      blocks <- lapply(starts,function(i){
        (i+1):(i+block.size)
      })
      if ( n_blocks != n_full_blocks ) {
        blocks <- c(blocks,list((n_full_blocks*block.size+1):num.epochs))
      }
      # randomize block order
      blocks <- blocks[sample.int(n_blocks)]
    }
    resampled_deltas <- unlist(lapply(blocks, function(indices){
      deltas[indices]
    }))
  } else {
    resampled_deltas <- sample(deltas,replace=replace)
  }
  
  x <- x0 + c(0,cumsum(resampled_deltas))
  
  if ( reverse.time ) {
    if ( keep.rate0 ) {
      warning("Cannot reverse time and keep first rate, setting reverse.time=FALSE")
    } else {
      x <- rev(x)
    }
  }
  
  if ( invert ) {
    x <- (x[1]^2) * 1/x
  }
  
  return (x)
}
