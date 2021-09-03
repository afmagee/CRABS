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
  l_derivative <- c(l_derivative[1], l_derivative)

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





pulled.speciation <- function( model ) {
  pulled.spec <- function(t, state, parameters){
    Lp <- state["Lp"]
    rp <- parameters["rp"]
    
    dLp = Lp * (model$p.delta(t) - Lp)
    return(list(dLp))
  }
  
  rho <- 1.0
  lambda0 <- model$lambda(0.0)
  
  parameters <- c(rp = model$p.delta)
  state <- c(Lp = rho*lambda0)
  
  res <- as.data.frame(radau(y = state, times = model$times, func = pulled.spec, parms = parameters), 
                        atol = 1e-06, rtol = 1e-06)
  
  # compute the derivatives
  #p_surv       <- p.survival.rate(v_spec0, v_ext0, times)
  #v_p_spec     <- v_spec0 * p_surv(times)
  Lp <- approxfun(res$time, res$Lp)
  
  return (Lp)
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

#' model2df
#'
#' @param model an object of class "ACDC"
#' @param gather boolean. Whether to return wide or long data frame
#'
#' @return a data frame
#' @export
#'
#' @examples
#' lambda <- function(t) 2.0 + sin(0.8*t)
#' mu <- function(t) 1.5 + exp(0.15*t)
#' times <- seq(from = 0, to = 4, length.out = 1000)
#' model <- congruence.class( lambda, mu, times = times)
#' 
#' model2df(model)
model2df <- function(model, gather = TRUE){
  times <- model$times
  l <- model$lambda(times)
  ex <- model$mu(times)
  
  df <- tibble::tibble("Time" = times,
               "Speciation" = l,
               "Extinction" = ex,
               "Net-diversification" = l - ex,
               "Relative extinction" = ex / l)
  if (gather){
    df <- tidyr::gather(df, "rate", "value", -Time)  
  }
  return(df)
}