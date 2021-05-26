R <- function(t, model, rtol){
  f <- function(u) model$lambda(u) - model$mu(u)
  quadgk(f, 0, t, tol = rtol)
}

logpsi <- function(s, t, model, rho, rtol){
  stopifnot(t > s)
  
  part_a <- R(t, model, rtol) - R(s, model, rtol)
  
  fb <- function(u) model$lambda(u)*exp(R(u, model, rtol)); fb <- Vectorize(fb)
  part_b <- log(1.0 + rho * quadgk(fb, 0.0, s, tol = rtol))
  
  fc <- function(u) model$lambda(u)*exp(R(u, model, rtol)); fc <- Vectorize(fc)
  part_c <- log(1.0 + rho * quadgk(fc, 0.0, t, tol = rtol))
  
  res = part_a + 2*(part_b - part_c)
}

fooE <- function(t, model, rtol){
  f <- function(s) model$lambda(s) * exp(R(s, model, rtol))
  f <- Vectorize(f)
  res <- exp(R(t, model, rtol)) / ((1.0 / rho) + quadgk(f, 0.0, t, tol = rtol))
}

#' Compute likelihood
#'
#' @param phy an object of class "phylo"
#' @param model an object of class "ACDC"
#' @param TESS whether or not to call `tess.likelihood()`
#' @param rtol relative tolerance for numerical integration
#' @param ... additional arguments passed to `tess.likelihood(...)`
#'
#' @return
#' @export
#'
#' @examples
#' library(ape)
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model <- congruence.class(lambda, mu, times = seq(0, 3, by = 0.005))
#' 
#' set.seed(123)
#' phy <- rcoal(25)
#' 
#' acdc.loglikelihood(phy, model)
acdc.loglikelihood <- function(phy, model, rho = 1.0, TESS = FALSE, rtol = 0.01, ...){
  if (TESS){
    th <- max(node.depth.edgelength(phy)) ## tree height
    
    ## translate the time coordinates
    mu <- function(t) model$mu(th - t)
    lambda <- function(t) model$lambda(th - t)
    times <- branching.times(phy)
    
    ## tess takes input as time starting at the root, increasing toward present
    res <- tess.likelihood(times = times, lambda = lambda, mu = mu, ...)
  }else{
    
    times <- sort(branching.times(phy), decreasing = TRUE)
    times <- unname(times)
    n <- length(times)
    
    res <- (n+1)*log(rho) + logpsi(0.0, times[1], model, rho, rtol) + log(model$lambda((times[1])))
    
    for (i in seq(from = 2, by = 1, to = n)){
      res <- res + log(model$lambda(times[i])) + logpsi(0.0, times[i], model, rho, rtol)
    }
    
    res <- res - log(fooE(times[1], model, rtol))
  }
  return(res)
}