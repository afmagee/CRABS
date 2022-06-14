congruent.speciation <- function( model, func.mu, ode_solver ) {
  
  lambda0 = model$lambda(0)
  times   = model$times
  v_p_div = sapply(times, model$p.delta )
  delta_t = model$delta_t
  
  v_ext1       <- sapply(times, func.mu)
  if(ode_solver){
    lambda_1   <- compute.speciation.ode(model, func.mu)
  }else{
    v_lambda1  <- compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
    lambda_1   <- approxfun( times, v_lambda1)
  }
  
  
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
  class(res) <- c("CRABS")


  
  return (res)
}


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



compute.speciation.ode <- function( model, mu ) {
  newLambda <- function(t, state, parameters){
    mu <- parameters["mu"]
    lambda <- parameters["lambda"]
    Lambda <- state["Lambda"]
    
    dLambda =  -Lambda^2 + Lambda *(model$p.delta(t) + mu(t))
    return(list(dLambda))
  }
  
  lambda0 <- model$lambda(0.0)
  
  parameters <- c(rp = model$p.delta,
                  mu = mu,
                  lambda = model$lambda)
  state <- c(Lambda = lambda0)
  
  res <- as.data.frame(deSolve::radau(y = state, times = model$times, func = newLambda, parms = parameters), 
                       atol = 1e-06, rtol = 1e-06)
  
  lambda <- approxfun(res$time, res$Lambda)
  
  return (lambda)
}
