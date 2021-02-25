###  create the pulled diversification rate
compute.pulled.diversification <- function( v_spec0, v_ext0, times, delta_t ) {

  # compute the derivatives
  l            <- v_spec0[-length(times)]
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


compute.extinction <- function( lambda0, v_p_div, v_ext1, delta_t ) {

  ### compute the new lambda
  v_lambda1    <- c()
  v_lambda1[1] <- lambda0

  for (j in 2:(NUM_TIME_DISCRETIZATIONS+1)) {

	tmp <- 4*v_lambda1[j-1]*delta_t + (v_p_div[j]*delta_t+v_ext1[j]*delta_t-1)^2
	v_lambda1[j] <- (sqrt(tmp) + v_p_div[j]*delta_t+v_ext1[j]*delta_t-1) / (2*delta_t)

  }

  return (v_lambda1)
}
