###  create the pulled diversification rate
compute.pulled.diversification <- function( v_spec0, v_ext0, delta_t ) {

  # compute the derivatives
  l            <- v_spec0[-length(v_ext0)]
  l_plus_one   <- v_spec0[-1]
  l_derivative <- (l_plus_one - l) / delta_t
  l_derivative <- c(l_derivative[1],l_derivative)

  # finally, add the 1/lambda * lambda dt to the pulled diversification rate
  v_p_div <- v_spec0 - v_ext0 + 1/v_spec0 * l_derivative

  return (v_p_div)
}



compute.speciation <- function( lambda0, v_p_div, v_ext1, delta_t ) {

  NUM_TIME_DISCRETIZATIONS = length(v_p_div)

  ### compute the new lambda
  v_lambda1    <- c()
  v_lambda1[1] <- lambda0

  for (j in 2:NUM_TIME_DISCRETIZATIONS) {

	tmp <- 4*v_lambda1[j-1]*delta_t + (v_p_div[j]*delta_t+v_ext1[j]*delta_t-1)^2
	v_lambda1[j] <- (sqrt(tmp) + v_p_div[j]*delta_t+v_ext1[j]*delta_t-1) / (2*delta_t)

  }

  return (v_lambda1)
}


compute.extinction <- function( v_p_div, v_spec1, delta_t ) {

  # compute the derivatives
  l            <- v_spec1[-length(v_spec1)]
  l_plus_one   <- v_spec1[-1]
  l_derivative <- (l_plus_one - l) / delta_t
  l_derivative <- c(l_derivative[1],l_derivative)

  # finally, add the 1/lambda * lambda dt to the pulled diversification rate
  v_mu1 <- v_spec1 - v_p_div + 1/v_spec1 * l_derivative

  return (v_mu1)
}


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
