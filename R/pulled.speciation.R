## From equation (57) in Louca & Pennell (2020) supplement. Page 10.
pulled.speciation <- function( model, rho = 1.0 ) {
  pulled.spec <- function(t, state, parameters){
    Lp <- state["Lp"]
    rp <- parameters["rp"]
    
    dLp = Lp * (model$p.delta(t) - Lp)
    return(list(dLp))
  }
  
  lambda0 <- model$lambda(0.0)
  
  parameters <- c(rp = model$p.delta)
  state <- c(Lp = rho*lambda0)
  
  res <- as.data.frame(deSolve::radau(y = state, times = model$times, func = pulled.spec, parms = parameters), 
                       atol = 1e-06, rtol = 1e-06)
  
  Lp <- approxfun(res$time, res$Lp)
  
  return (Lp)
}
