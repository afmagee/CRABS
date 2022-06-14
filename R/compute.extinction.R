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
  
  v_spec1        <- sapply(times, func.lambda)
  #v_mu1          <- compute.extinction( v_p_div, v_spec1, delta_t )
  lambda_1       <- func.lambda
  v_mu1          <- compute.extinction( lambda_1, v_spec1, v_p_div, times )
  mu_1           <- approxfun( times, v_mu1)
  
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
  class(res) <- c("CRABS")
  
  return (res)
}

compute.extinction <- function(lambda_1, v_spec1, v_p_div, times) {
  back <- pracma::fderiv(lambda_1, times, method = "backward")
  forw <- pracma::fderiv(lambda_1, times, method = "forward")
  m <- rbind(forw, back)
  lambda_deriv <- apply(m, 2, mean, na.rm = TRUE)
  
  v_mu1 <- v_spec1 - v_p_div + 1/v_spec1 * lambda_deriv
}
# 
# compute.extinction <- function( v_p_div, v_spec1, delta_t ) {
# 
#   # compute the derivatives
#   l            <- v_spec1[-length(v_spec1)]
#   l_plus_one   <- v_spec1[-1]
#   l_derivative <- (l_plus_one - l) / delta_t
#   l_derivative <- c(l_derivative[1],l_derivative)
# 
#   # finally, add the 1/lambda * lambda dt to the pulled diversification rate
#   v_mu1 <- v_spec1 - v_p_div + 1/v_spec1 * l_derivative
# 
#   return (v_mu1)
# }
