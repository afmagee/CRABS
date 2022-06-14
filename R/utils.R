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
#' @param model an object of class "CRABS"
#' @param gather boolean. Whether to return wide or long data frame
#' @param rho the sampling fraction at the present. Used to calculate the pulled speciation rate
#' @param compute.pulled.rates whether to compute the pulled rates
#'
#' @return a data frame
#' @export
#'
#' @examples
#' lambda <- function(t) 2.0 + sin(0.8*t)
#' mu <- function(t) 1.5 + exp(0.15*t)
#' times <- seq(from = 0, to = 4, length.out = 1000)
#' model <- create.model( lambda, mu, times = times)
#' 
#' model2df(model)
model2df <- function(model, gather = TRUE, rho = 1.0, compute.pulled.rates = TRUE){
  times <- model$times
  l <- model$lambda(times)
  ex <- model$mu(times)
  
  if (compute.pulled.rates){
    p.lambda <- pulled.speciation(model, rho = rho)(times)
    p.delta <- model$p.delta(times)  
    p.mu <- model$lambda(0.0) - p.delta
  }else{
    p.lambda <- NULL
    p.delta <- NULL
    p.mu <- NULL
  }
  
  
  df <- tibble::tibble("Time" = times,
               "Speciation" = l,
               "Extinction" = ex,
               "Net-diversification" = l - ex,
               "Relative extinction" = ex / l,
               "Pulled net-diversification" = p.delta,
               "Pulled speciation" = p.lambda,
               "Pulled extinction" = p.mu)
  if (gather){
    df <- gather(df, "rate", "value", -Time)  
  }
  return(df)
}