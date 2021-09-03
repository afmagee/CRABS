#' Create the piecewise-constant speciation rate from extinction rate
#'
#' @param congruence.class The model congruence class
#' @param func.mu A extinction rate function
#' @param list.funcs.mu A list of extinction rate functions
#' @return An object with a list of speciation and extinction rate functions.
#' @export
ACDC.compute.speciation <- function( congruence.class, func.mu=NULL, list.funcs.mu=NULL ) {

  lambda0 = congruence.class$lambda(0)
  times   = congruence.class$times
  v_p_div = sapply(times, congruence.class$p.delta )
  delta_t = congruence.class$delta_t

  lambda_1 = list()
  mu_1     = list()
  if ( is.null( func.mu ) == FALSE ) {
     v_ext1         = sapply(times, func.mu)#func.mu( times )
     v_lambda1      = compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
     lambda_1[[1]]  = approxfun( times, v_lambda1)
     mu_1[[1]]      = func.mu
  } else if ( is.null( list.funcs.mu ) == FALSE  ) {
    for (i in 1:length(list.funcs.mu)) {
       v_ext1         = list.funcs.mu[[i]]( times )
       v_lambda1      = compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
       lambda_1[[i]]  = approxfun( times, v_lambda1)
       mu_1[[i]]      = list.funcs.mu[[i]]
    }
  }


  res = list( lambda0=congruence.class$lambda,
              mu0=congruence.class$mu,
              lambda1=lambda_1,
              mu1=mu_1,
              max.t=congruence.class$max.t )

  return (res)
}

#' Create the congruent model from extinction rate
#'
#' @param congruence.class The model congruence class
#' @param func.mu A extinction rate function
#' @param list.funcs.mu A list of extinction rate functions
#' @return An object with a list of speciation and extinction rate functions.
#' @export
congruent.speciation <- function( model, func.mu ) {
  
  lambda0 = model$lambda(0)
  times   = model$times
  v_p_div = sapply(times, model$p.delta )
  delta_t = model$delta_t
  
  v_ext1         = sapply(times, func.mu)#func.mu( times )
  v_lambda1      = compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
  lambda_1  = approxfun( times, v_lambda1)
  mu_1      = func.mu
  if (length(mu_1(times)) == 1){
    mu_1 = Vectorize(mu_1)
  }
  
  ## create the parameter transformations as rate functions
  func_div    <- function(t) lambda_1(t) - mu_1(t)
  func_turn   <- function(t) mu_1(t) / lambda_1(t)
  
  #stop()
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

  #stop()
  for (j in 2:NUM_TIME_DISCRETIZATIONS) {
  # Finite forwards difference
    v_lambda1[j] <- v_lambda1[j-1]*(1 + delta_t * (v_p_div[j-1] - v_lambda1[j-1] + v_ext1[j-1]))
    
  # Finite backward difference
	#tmp <- 4*v_lambda1[j-1]*delta_t + (v_p_div[j]*delta_t+v_ext1[j]*delta_t-1)^2
	#v_lambda1[j] <- (sqrt(tmp) + v_p_div[j]*delta_t+v_ext1[j]*delta_t-1) / (2*delta_t)

  }
  #stop()

  return (v_lambda1)
}
