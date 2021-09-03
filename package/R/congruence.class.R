#' Computes the congruent class, i.e., the pulled rates.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param num.intervals The number of discretization.
#' @return A list of rate functions representing this congruence class.
#' @export
#' @examples
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model <- congruence.class(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' model
congruence.class <- function(func_spec0, func_ext0, times = seq(from = 0, to = 5, by = 0.005), func_p_spec = NULL, func_p_div = NULL) {

  ## create our vector of times (i.e., change-points)
  ## for the piecewise linear approximation
  #times       <- (0:num.intervals) / num.intervals * max.t
  
  max.t <- max(times)
  num.intervals <- length(times)
  delta_t <- times[2] - times[1]

  ## create the vector of rate values at the change points
  v_spec0     <- func_spec0( times )
  if(length(v_spec0) == 1){
    v_spec0 <- rep(v_spec0, num.intervals)
    func_spec0 <- Vectorize(func_spec0)
  }
  v_ext0      <- func_ext0( times )
  if(length(v_ext0) == 1){
    v_ext0 <- rep(v_ext0, num.intervals)
    func_ext0 <- Vectorize(func_ext0)
  }

  ## create the parameter transformations as rate functions
  func_div    <- function(t) func_spec0(t) - func_ext0(t)
  func_turn   <- function(t) func_ext0(t) / func_spec0(t)
  

  ## compute the pulled diversification rate
  if(missing("func_p_div")){
    v_p_div     <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )
    func_p_div  <- approxfun(times,v_p_div)
  }
  
  ## compute the pulled speciation rate
  # if(missing("func_p_spec")){
  #   print(microbenchmark(v_p_spec    <- compute.pulled.speciation( v_spec0, v_ext0, times )))
  #   v_p_spec    <- compute.pulled.speciation( v_spec0, v_ext0, times )
  #   func_p_spec <- approxfun(times,v_p_spec)
  # }
  
  res = list(lambda=func_spec0,
             mu=func_ext0,
             delta=func_div,
             epsilon=func_turn,
             p.delta=func_p_div,
             #p.lambda=func_p_spec,
             times=times,
             max.t = max.t,
             delta_t = delta_t,
             num.intervals = num.intervals)
  class(res) <- c("ACDC")

  return (res)
}

#' Computes the congruent class, i.e., the pulled rates.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param num.intervals The number of discretization.
#' @return A list of rate functions representing this congruence class.
#' @export
#' @examples
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model <- congruence.class(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' model
create.model <- function(func_spec0, func_ext0, times = seq(from = 0, to = 5, by = 0.005), func_p_spec = NULL, func_p_div = NULL) {
  
  ## create our vector of times (i.e., change-points)
  ## for the piecewise linear approximation
  max.t <- max(times)
  num.intervals <- length(times)
  delta_t <- times[2] - times[1]
  
  ## create the vector of rate values at the change points
  v_spec0     <- func_spec0( times )
  if(length(v_spec0) == 1){
    v_spec0 <- rep(v_spec0, num.intervals)
    func_spec0 <- Vectorize(func_spec0)
  }
  v_ext0      <- func_ext0( times )
  if(length(v_ext0) == 1){
    v_ext0 <- rep(v_ext0, num.intervals)
    func_ext0 <- Vectorize(func_ext0)
  }
  
  ## create the parameter transformations as rate functions
  func_div    <- function(t) func_spec0(t) - func_ext0(t)
  func_turn   <- function(t) func_ext0(t) / func_spec0(t)
  
  ## compute the pulled diversification rate
  if(missing("func_p_div")){
    v_p_div     <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )
    func_p_div  <- approxfun(times,v_p_div)
  }
  
  res = list(lambda=func_spec0,
             mu=func_ext0,
             delta=func_div,
             epsilon=func_turn,
             p.delta=func_p_div,
             times=times,
             max.t = max.t,
             delta_t = delta_t,
             num.intervals = num.intervals)
  class(res) <- c("ACDC")
  
  return (res)
}


#' Print method for ACDC object
#'
#' @param model and object of class ACDC
#' @param ... other arguments
#'
#' @return
#'
#' @examples
#' @export
print.ACDC <- function(model, ...){
  cat("Piecewise-linear birth-death model\n")
  cat("Knots:", length(model$times), "\n")
  cat("Delta-tau:", model$delta_t, "\n")
  ACDC.plot.rates(model)
  invisible()
}