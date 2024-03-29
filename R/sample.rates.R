#' Sample custom functions through time.
#'
#' @param times the time knots
#' @param lambda0 The rate at present
#' @param rsample Function to sample next rate
#' @param rsample0 Function to sample rate at present
#' @param autocorrelated Should rates be autocorrelated?
#' @return Sampled rate vector
#' @export
#' @examples
#' data("primates_ebd")
#' 
#' l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' times <- primates_ebd[["time"]]
#' 
#' model <- create.model(l, mu, times)
#' 
#' rsample <- function(n) runif(n, min = 0.0, max = 0.9)
#' mu <- sample.rates(times, 0.5, rsample = rsample)
#' 
#' 
#' model_set <- congruent.models(model, mus = mu)
#' 
#' model_set
sample.rates <- function(times, 
                         lambda0=NULL, 
                         rsample=NULL, 
                         rsample0=NULL,
                         autocorrelated=FALSE) {
  num.epochs <- length(times)
  
  N_SAMPLES = ifelse( is.null(lambda0), num.epochs, num.epochs -1)
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
  func_rates <- approxfun(times, new_rates)
  
  return (func_rates)
}


#' Samples simple increase/decrease models through time with noise.
#'
#' @param times the time knots
#' @param rate0 The rate at present, otherwise drawn randomly.
#' @param model "MRF" for pure MRF model, otherwise MRF has a trend of type "exponential","linear", or "episodic<n>"
#' @param direction "increase" or "decrease" (measured in past to present)
#' @param noisy If FALSE, no MRF noise is added to the trajectory
#' @param MRF.type "HSMRF" or "GMRF", type for stochastic noise.
#' @param monotonic Whether the curve should be forced to always move in one direction.
#' @param fc.mean Determines the average amount of change when drawing from the model.
#' @param rate0.median When not specified, rate at present is drawn from a lognormal distribution with this median.
#' @param rate0.logsd When not specified, rate at present is drawn from a lognormal distribution with this sd
#' @param mrf.sd.scale scale the sd of the mrf process up or down. defaults to 1.0
#' @param min.rate The minimum rate (rescaling fone after after drawing rates).
#' @param max.rate The maximum rate (rescaling fone after after drawing rates).
#' @return Speciation or extinction rate at a number of timepoints.
#' @export
#' @examples
#' data("primates_ebd")
#' 
#' l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' times <- primates_ebd[["time"]]
#' 
#' model <- create.model(l, mu, times)
#' 
#' mus <- sample.basic.models(times = times, 
#'                                rate0 = 0.05, 
#'                                "MRF", 
#'                                MRF.type = "HSMRF", 
#'                                fc.mean = 2.0, 
#'                                min.rate = 0.0, 
#'                                max.rate = 1.0)
#' 
#' model_set <- congruent.models(model, mus = mus)
#' 
#' model_set
sample.basic.models <- function(times, 
                                rate0=NULL, 
                                model="exponential", 
                                direction="decrease", 
                                noisy=TRUE, 
                                MRF.type="HSMRF", 
                                monotonic=FALSE, 
                                fc.mean=3, 
                                rate0.median=0.1, 
                                rate0.logsd=1.17481, 
                                mrf.sd.scale = 1.0,
                                min.rate=0, 
                                max.rate=10) {
  num.epochs <- length(times)
  # recover()
  
  # We use rejection sampling to find a model that fits within minimum and maximum rates
  # To speed up sampling, we break this into three rejection sampling steps, handling separately the:
  #    1) rate at present
  #    2) fold change (given rate at present)
  #    3) MRF noise
  
  # rate at present
  x0 <- min.rate - 10
  if ( !is.null(rate0) ) {
    if ( rate0 < min.rate || rate0 > max.rate ) {
      stop("User-defined rate0 is outside [min.rate,max.rate].")
    }
    x0 <- rate0
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
  
  if ( direction == "increase" ) {
    fc <- 1/fc
  } else if ( direction != "decrease" ) {
    stop("Invalid \"direction\"")
  }
  
  # Deterministic component of trajectory
  num_deltas <- num.epochs -1
  x <- numeric(num.epochs)
  x[1] <- x0
  if ( model == "exponential" ) {
    delta_deterministic <- rep(log(fc)/(num_deltas),num_deltas)
    x[2:(num.epochs)] <- x[1] * exp(cumsum(delta_deterministic))
  } else if ( model == "linear" ) {
    delta_deterministic <- rep(((x0*fc)-x0)/(num_deltas),num_deltas)
    x[2:(num.epochs)] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("episodic",model) ) {
    delta_deterministic <- rep(0,num_deltas)
    njumps <- as.numeric(gsub("episodic","",model)) - 1
    if (njumps < 1) {
      stop("Too few episodes in episodic model")
    }
    if ( njumps == 1 ) {
      delta_deterministic[sample.int(num_deltas,1)] <- (fc * x0) - x0
    } else {
      delta_deterministic[sample.int(num_deltas,njumps)] <- ((fc * x0) - x0) * rdirichlet(njumps,1)
    }
    x[2:(num.epochs)] <- x[1] + cumsum(delta_deterministic)
  } else if ( grepl("MRF",model) ) {
    x <- rep(x0,num.epochs)
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
      delta_stochastic <- rep(Inf,num_deltas)
      noise <- rep(Inf,num_deltas)
      # Avoid infinities, this loop should rarely trigger
      while ( any(!is.finite(noise)) ) {
        if ( MRF.type == "HSMRF"  || MRF.type == "HSRF") {
          zeta <- get.hsmrf.global.scale(num.epochs)
          gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
          sigma <- abs(rcauchy(num_deltas,0,1))
          delta_stochastic <- rnorm(num_deltas,0,sigma*gamma*zeta*mrf.sd.scale)
        } else if ( MRF.type == "GMRF" ) {
          zeta <- get.gmrf.global.scale(num.epochs)
          gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
          delta_stochastic <- rnorm(num_deltas,0,gamma*zeta*mrf.sd.scale)
        } else {
          stop("Invalid \"MRF.type\"")
        }
        
        if ( monotonic ) {
          delta_stochastic <- abs(delta_stochastic)
          
          if ( direction == "increase" ) {
            delta_stochastic <- -delta_stochastic
          }
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
  
  func_rates <- approxfun(times, rates)
  
  return (func_rates)
}


#' Jointly samples speciation and extinction trajectories through time, with noise.
#'
#' @param times the time knots
#' @param p.delta The The pulled diversification rate function (measured in time before present).
#' @param lambda0 The speciation rate at present.
#' @param mu0 The extinction rate at present, otherwise drawn randomly.
#' @param MRF.type "HSMRF" or "GMRF", type for stochastic noise.
#' @param beta.param Parameters of the Beta distribution used for 
#' @param mu0.median When not specified, extinction rate at present is drawn from a lognormal distribution with this median.
#' @param mu0.logsd When not specified, extinction rate at present is drawn from a lognormal distribution with this sd
#' @param mrf.sd.scale scale the sd of the mrf process up or down. defaults to 1.0
#' @param min.lambda The minimum speciation rate (rescaling done after after drawing rates).
#' @param min.mu The minimum extinction rate (rescaling done after after drawing rates).
#' @param max.lambda The maximum speciation rate (rescaling done after after drawing rates).
#' @param max.mu The maximum extinction rate (rescaling done after after drawing rates).
#' @param min.p The lower bound of parameter p's trajectory.
#' @param max.p The upper bound of parameter p's trajectory.
#' @return Speciation or extinction rate at a number of timepoints.
#' @export
#' @examples
#' data("primates_ebd")
#' 
#' l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' times <- primates_ebd[["time"]]
#' 
#' model <- create.model(l, mu, times)
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
sample.basic.models.joint <- function(times, 
                                      p.delta,
                                      lambda0,
                                      mu0=NULL, 
                                      MRF.type="HSMRF", 
                                      beta.param=c(0.3,0.3),
                                      mu0.median=0.1, 
                                      mu0.logsd=1.17481, 
                                      mrf.sd.scale = 1.0,
                                      min.lambda=0, 
                                      min.mu=0, 
                                      max.lambda=10, 
                                      max.mu=10,
                                      min.p=-0.05, 
                                      max.p=1.05) {
  num.epochs <- length(times)
  delta_t <- times[2]-times[1]

  # rate at present
  if ( lambda0 < min.lambda || lambda0 > max.lambda ) {
    stop("User-defined lambda0 is outside [min.lambda,max.lambda].")
  }
  x0 <- lambda0
  
  y0 <- min.mu - 10
  if ( !is.null(mu0) ) {
    if ( mu0 < min.mu || mu0 > max.mu ) {
      stop("User-defined mu0 is outside [min.mu,max.mu].")
    }
    y0 <- mu0
  } else {
    while ( y0 < min.mu || y0 > max.mu ) {
      y0 <- rlnorm(1,log(mu0.median),mu0.logsd)
    }
  }
  
  lambdas <- rep(x0, num.epochs)
  mus <- rep(y0, num.epochs)
  
  ### MRF joint trajectories between λ_{i-1} and λ* (:= λ such as μ_i=μ_{i-1}) ###
    
  in_bounds = FALSE
  while (!in_bounds){
    in_bounds = TRUE
    
    # Generate the MRF trajectory
    if (MRF.type == "HSMRF"  || MRF.type == "HSRF") {
      zeta <- get.hsmrf.global.scale(num.epochs)
      gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
      sigma <- abs(rcauchy(num.epochs-1,0,1))
      delta_stochastic <- rnorm(num.epochs-1,0,sigma*gamma*zeta*mrf.sd.scale)
    } else if (MRF.type == "GMRF") {
      zeta <- get.gmrf.global.scale(num.epochs)
      gamma <- min(abs(rcauchy(1,0,1)),1000) # avoid numerical instability
      delta_stochastic <- rnorm(num.epochs-1,0,gamma*zeta*mrf.sd.scale)
    } else {
      stop("Invalid \"MRF.type\"")
    }
    trajectory <- cumsum(delta_stochastic)
    
    # Fix a random point in the trajectory at a value drawn from a Beta distribution
    trajectory <- trajectory + rbeta(1, beta.param[1], beta.param[2]) - sample(trajectory, 1)
    # Rescale or slide the trajectory if it exceeds the allowed bounds
    trajectory <- trajectory / max((max(trajectory)-min(trajectory))/(max.p-min.p), 1)
    trajectory <- trajectory + min.p - min(min(trajectory), min.p)
    trajectory <- trajectory + max.p - max(max(trajectory), max.p)
    trajectory <- c(NA, trajectory)
    
    for (i in 2:num.epochs){
      # compute the λ_min and λ* parameters
      lambda_min <- (p.delta(times[i])-1/delta_t + sqrt((1/delta_t-p.delta(times[i]))**2+4*lambdas[i-1]/delta_t))/2
      b = p.delta(times[i]) + mus[i-1] - 1/delta_t
      lambda_star = (b + sqrt(b**2+4*lambdas[i-1]/delta_t))/2
      # get λ_i between  λ_{i-1} and λ* according to the trajectory
      lambda_i <- lambdas[i-1] + trajectory[i]*(lambda_star-lambdas[i-1])
      # force λ_i to be higher than  λ_min
      if (lambda_i<lambda_min){
        lambda_i <- lambda_min
        mu_i <- 0
      }else{
        mu_i <- lambda_i - p.delta(times[i]) + (lambda_i-lambdas[i-1])/(lambda_i*delta_t)
      }
      
      if (mu_i > max.mu || lambda_i > max.lambda){
        in_bounds = FALSE
        break
      }
      lambdas[i] <- lambda_i
      mus[i] <- mu_i
    }
  }

  func_lambdas <- approxfun(times, lambdas)
  func_mus <- approxfun(times, mus)
  
  return (list(func_lambdas=func_lambdas, func_mus=func_mus))
}
