#' Create the congruent model from extinction rate
#'
#' @param model The reference model
#' @param func.mu An extinction rate function
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
  
  res <- list(lambda = lambda_1,
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

  for (j in 2:NUM_TIME_DISCRETIZATIONS) {
    # Finite backward difference
  	tmp <- 4*v_lambda1[j-1]*delta_t + (v_p_div[j]*delta_t+v_ext1[j]*delta_t-1)^2
  	v_lambda1[j] <- (sqrt(tmp) + v_p_div[j]*delta_t+v_ext1[j]*delta_t-1) / (2*delta_t)
  }

  return (v_lambda1)
}
