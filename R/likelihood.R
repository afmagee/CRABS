dLTT <- function(model, times, N0, rho = 1.0){
  lambda0 <- model$lambda(0.0)
  
  odeN <- function(time, state, parameters, lambda, mu){
    l <- lambda(time)
    m <- mu(time)
    with(as.list(c(state, parameters)), {
      dN = N * (m - l)
      dE = m - E * (l + m) + E^2 * l
      return(list(c(dN, dE), l, m))
    })
  }
  
  parameters <- list()
  state <- c("N" = N0 / rho, "E" = 1 - rho)
  
  res <- as.data.frame(deSolve::radau(y = state, times = times, 
                                      func = odeN, parms = parameters, 
                                      lambda = model$lambda, mu = model$mu))
  M <- res$N * (1 - res$E)
  return(M)
}

#' Compute likelihood
#'
#' @param phy an object of class "phylo"
#' @param model an object of class "CRABS"
#' @param rho the taxon sampling fraction
#'
#' @return the log-likelihood of the tree given the model
#' @export
#'
#' @examples
#' library(ape)
#' lambda <- function(t) exp(0.3*t) - 0.5*t
#' mu <- function(t) exp(0.3*t) - 0.2*t - 0.8
#'  
#' model <- create.model(lambda, mu, times = seq(0, 3, by = 0.005))
#' 
#' set.seed(123)
#' phy <- rcoal(25)
#' 
#' crabs.loglikelihood(phy, model)
crabs.loglikelihood <- function(phy, model, rho = 1.0){
  times <- model$times
  N0 <- length(phy$tip.label)
  M <- approxfun(times, dLTT(model, times, N0, rho))
  
  bt <- branching.times(phy)
  n <- length(bt)
  
  dM = tail(fderiv(M, bt), n = -1)
  
  logL <- log(M(bt[1])) - (n+1)*log(M(0.0)) + sum(log(-dM))
  return(logL)
}