#' Computes the congruent class, i.e., the pulled rates.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param num.intervals The number of discretization.
#' @return A list of rate functions representing this congruence class.
#' @export
#' @examples
#' #TODO
congruence.class <- function(func_spec0, func_ext0, max.t, num.intervals=1000 ) {

  ## create our vector of times (i.e., change-points)
  ## for the piecewise linear approximation
  times     <- (0:num.intervals) / num.intervals * max.t
  delta_t   <- max.t / num.intervals

  ## create the vector of rate values at the change points
  v_spec0   <- func_spec0( times )
  v_ext0    <- func_ext0( times )

  ## create the parameter transformations as rate functions
  func_div  <- function(t) func_spec0(t) - func_ext0(t)
  func_turn <- function(t) func_ext0(t) / func_spec0(t)

  ## compute the pulled diversification rate
  v_p_div    <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )
  func_p_div <- approxfun(times,v_p_div)

  ## compute the pulled speciation rate
  v_p_spec   <- compute.pulled.speciation( v_spec0, v_ext0, times )
  func_p_spec <- approxfun(times,v_p_spec)

  res = list(lambda=func_spec0,
             mu=func_ext0,
             delta=func_div,
             epsilon=func_turn,
             p.delta=func_p_div,
             p.lambda=func_p_spec,
             max.t=max.t)

  return (res)
}
