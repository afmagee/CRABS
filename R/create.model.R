#' Computes the congruent class, i.e., the pulled rates.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param times the time knots for the piecewise-linear rate functions
#' @param func_p_spec the pulled speciation rate function
#' @param func_p_div the pulled net-diversification rate function
#' @return A list of rate functions representing this congruence class.
#' @export
#' @examples
#' lambda1 <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu1 <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model1 <- create.model(lambda1, mu1, times = seq(0, 5, by = 0.005))
#' 
#' model1
#' 
#' data("primates_ebd")
#' 
#' lambda2 <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu2 <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' model2 <- create.model(lambda2, mu2, primates_ebd[["time"]])
#' 
#' model2
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
    #v_p_div     <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )
    back <- pracma::fderiv(func_spec0, times, method = "backward")
    forw <- pracma::fderiv(func_spec0, times, method = "forward")
    m <- rbind(forw, back)
    lambda_deriv <- apply(m, 2, mean, na.rm = TRUE)
    
    v_p_div <- v_spec0 - v_ext0 + (1/v_spec0) * lambda_deriv
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
  class(res) <- c("CRABS")
  
  return (res)
}