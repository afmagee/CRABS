## This carries out integration for the time-varying
## speciation/extinction model with functions lambda and mu,
## outputting at the times in the vector 'times'.
p.survival.rate <- function(v_lambda, v_mu, times) {

  lambda   <- approxfun(times,v_lambda)
  mu       <- approxfun(times,v_mu)

  E <- array(0,length(times))

  dt           <- (times[2] - times[1]) / 100
  current_time <- 0
  this_E <- 0
  for (i in 1:length(times)) {
    while ( current_time < times[i] ) {
      d             <- mu(current_time)
      b             <- lambda(current_time)
      this_E        <- this_E + ( d - (b+d)*this_E + b*this_E*this_E)*dt
      current_time  <- current_time + dt
    }
    E[i] <- this_E
  }

  p_surv <- approxfun( times, 1-E )

  return(p_surv)
}
