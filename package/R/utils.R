#' Create the piecewise-constant pulled diversification rate
#'
#' @param v_spec0 The speciation rate at all changepioints
#' @param v_ext0 The extinction rate at all changepoints
#' @param delta_t The width of each grid cell
#' @return Pulled diversification rate at all changepoints
#' @keywords internal
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


#' Create the piecewise-constant speciation rate
#'
#' @param lambda0 The rate at present
#' @param v_p_div The pulled diversification rate at all changepioints
#' @param v_ext1 The extinction rate at all changepoints
#' @param delta_t The width of each grid cell
#' @return Speciation rate at all changepoints
#' @keywords internal
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

#' Create the piecewise-constant extinction rate
#'
#' @param lambda0 The rate at present
#' @param v_p_div The pulled diversification rate at all changepioints
#' @param v_ext1 The extinction rate at all changepoints
#' @param delta_t The width of each grid cell
#' @return Extinction rate at all changepoints
#' @keywords internal
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

#' One draw from a dirichlet distribution with concentration parameter alpha.
#'
#' @param dim The dimension
#' @param alpha Concentration parameter
#' @return Vector
#' @keywords internal
rdirichlet <- function(dim,alpha) {
  x <- rgamma(dim,alpha,1)
  return(x/sum(x))
}
