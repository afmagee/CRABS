## This carries out integration for the time-varying
## speciation/extinction model with functions lambda and mu,
## outputting at the times in the vector 'times'.
p.survival.rate <- function(lambda, mu, times) {
  obj <- function(t, state, parms) {
    r <- mu(t) - lambda(t)
    list(c(r))
  }

  tmax <- times[length(times)]
  t.crit <- c(0,tmax)
  n <- 1
  bin <- findInterval(times, t.crit, TRUE)
  y <- c(0)

  xx <- c()
  r_points <- c()
  # we compute the integrals stepwise from t.crit[i] to t.crit[i+1]
  for ( i in seq_len(n) ) {
    j <- i + 1L
    ti <- c(t.crit[[i]], times[bin == i], t.crit[[j]]-1E-9)
    yi <- lsoda(y, ti, obj, tcrit=t.crit[[j]])

    xx <- c(xx,ti)
    r_points <- c(r_points,yi[,2])

    y <- yi[nrow(yi),-1]

  }

  # add a final point slightly over the given time interval. we need this so that some of the numerical routines do not break when calling r(t).
  xx <- c(-1E-4 ,xx, tmax + 1e-8)
  r_points <- c(0,r_points,r_points[length(r_points)])
  rate <- approxfun(xx,r_points)

  # compute the integral int_{x}^{T} mu(t)*exp(r(t)) dt
  f <- function(t) pmin(mu(t)*exp(rate(t)),1E100)
  probs <- array(0,length(xx))
  for (i in length(xx):2) {
    u <- xx[i]
    l <- xx[i-1]
    val <- integrate(f,upper=u,lower=l)$value
    probs[i-1] <- val + probs[i]

  }
  surv <- approxfun(xx,probs)

  list(r=rate,s=surv)

  p_surv <- function( t ) {
    den <- ( 1 + ( ( surv(0) - surv(t) ) / exp(rate(0)) ) )
    return ( 1/den )
  }

  return(p_surv)
}
