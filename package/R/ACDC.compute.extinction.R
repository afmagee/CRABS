#' Create the piecewise-constant extinction rate from speciation rate
#'
#' @param congruence.class The model congruence class
#' @param func.lambda A speciation rate function
#' @param list.funcs.lambda A list of speciation rate functions
#' @return An object with a list of speciation and extinction rate functions.
#' @export
ACDC.compute.extinction <- function( congruence.class, func.lambda=NULL, list.funcs.lambda=NULL ) {

  times   = congruence.class$times
  v_p_div = sapply(times, congruence.class$p.delta)#congruence.class$p.delta( times )
  delta_t = congruence.class$delta_t

  lambda_1 = list()
  mu_1     = list()
  #stop()
  if ( is.null( func.lambda ) == FALSE ) {
     ## make sure that the initial conditions holds
     if ( abs(congruence.class$lambda(0) - func.lambda(0)) > 1E-8 ) {
       stop("The initial values of the reference and alternative speciation rate functions are not identical")
     }
     v_spec1        = sapply(times, func.lambda)
     v_mu1          = compute.extinction( v_p_div, v_spec1, delta_t )
     lambda_1[[1]]  = func.lambda
     mu_1[[1]]      = approxfun( times, v_mu1)
  } else if ( is.null( list.funcs.lambda ) == FALSE  ) {
    for (i in 1:length(list.funcs.lambda)) {
       ## make sure that the initial conditions holds
       if ( abs(congruence.class$lambda(0) - list.funcs.lambda[[i]](0)) > 1E-8 ) {
         stop("The initial values of the reference and alternative speciation rate functions are not identical")
       }
       v_spec1        = list.funcs.lambda[[i]]( times )
       v_mu1          = compute.extinction( v_p_div, v_spec1, delta_t )
       lambda_1[[i]]  = list.funcs.lambda[[i]]
       mu_1[[i]]      = approxfun( times, v_mu1)
    }
  }


  res = list( lambda0=congruence.class$lambda,
              mu0=congruence.class$mu,
              lambda1=lambda_1,
              mu1=mu_1,
              max.t=congruence.class$max.t )

  return (res)
}

#' @export
congruent.extinction <- function( model, func.lambda ) {
  times   = model$times
  v_p_div = sapply(times, model$p.delta)
  delta_t = model$delta_t
  
  if (length(func.lambda(times)) == 1){
    func.lambda <- Vectorize(func.lambda)
  }
  
  ## make sure that the initial conditions holds
  if ( abs(model$lambda(0) - func.lambda(0)) > 1E-8 ) {
    stop("The initial values of the reference and alternative speciation rate functions are not identical")
  }
  
  v_spec1        = sapply(times, func.lambda)
  v_mu1          = compute.extinction( v_p_div, v_spec1, delta_t )
  lambda_1  = func.lambda
  mu_1      = approxfun( times, v_mu1)
  
  ## create the parameter transformations as rate functions
  func_div    <- function(t) lambda_1(t) - mu_1(t)
  func_turn   <- function(t) mu_1(t) / lambda_1(t)
  
  res = list(lambda = lambda_1,
             mu = mu_1,
             delta = func_div,
             epsilon = func_turn,
             p.delta = model$p.delta,
             times = times,
             max.t = model$max.t,
             delta_t = model$delta_t,
             num.intervals = model$num.intervals)
  class(res) <- c("ACDC")
  
  return (res)
}

#' Create the piecewise-constant extinction rate
#'
#' @param lambda0 The rate at present
#' @param v_p_div The pulled diversification rate at all change-points
#' @param v_ext1 The extinction rate at all change-points
#' @param delta_t The width of each grid cell
#' @return Extinction rate at all change-points
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
