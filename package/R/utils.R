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


#' Create the piecewise-constant pulled speciation rate
#'
#' @param v_spec0 The speciation rate at all change-points
#' @param v_ext0 The extinction rate at all change-points
#' @param times The change-points
#' @return Pulled speciation rate at all change-points
#' @keywords internal
compute.pulled.speciation <- function( v_spec0, v_ext0, times ) {

  # compute the derivatives
  p_surv       <- p.survival.rate(v_spec0, v_ext0, times)
  v_p_spec     <- v_spec0 * p_surv(times)

  return (v_p_spec)
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

#' Global scale for HSMRF using linear interpolation of pre-computed values
#'
#' @param v The number of pieces in the approximation
#' @return Global scale
#' @keywords internal
get.hsmrf.global.scale <- approxfun(x=c(2,10,20,50,100,200,500,1000,2000,5000,10000,100000),
                                    y=c(3.971988,0.04446432,0.01653494,0.005024999,0.002093737,0.0008441697,0.0002137604,4.947984e-05,4.518384e-10,7.099068e-11,3.113958e-11,2.941689e-12))


#' Global scale for GMRF using linear interpolation of pre-computed values
#'
#' @param v The number of pieces in the approximation
#' @return Global scale
#' @keywords internal
get.gmrf.global.scale <- approxfun(x=c(2,10,20,50,100,200,500,1000,2000,5000,10000,100000),
                                    y=c(2.871193,0.1064935,0.04975563,0.01911503,0.009376335,0.004549474,0.001693932,0.0007408964,0.0002640923,0.0001002221,7.352401e-05,4.42448e-05))
