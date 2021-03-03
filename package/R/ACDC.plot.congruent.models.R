#' Plots the rate functions including the pulled rates.
#'
#' @param obj A list of rate functions produced by .
#' @return void.
#' @export
#' @examples
#' #TODO
ACDC.plot.congruent.models <- function( obj ) {

  ## general settings
  num.intervals = 10000
  times         = (0:num.intervals) / num.intervals * obj$max.t
  num.models    = length( obj$lambda1 )

  ## plot settings
  col_lambda    = brewer.pal(num.models+1, "Blues")[2:(num.models+1)]
  col_mu        = brewer.pal(num.models+1, "Reds")[2:(num.models+1)]
  col_delta     = brewer.pal(num.models+1, "Purples")[2:(num.models+1)]
  col_epsilon   = brewer.pal(num.models+1, "Greens")[2:(num.models+1)]
  this.lwd      = 1

  par(mfrow=c(2,2))

  lambda = obj$lambda0
  Y_MIN <- min( lambda(times) )
  Y_MAX <- max( lambda(times) )
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    Y_MIN <- min(Y_MIN, lambda(times))
    Y_MAX <- max(Y_MAX, lambda(times))
  }
  plot( NA, NA, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    curve(lambda, lwd=this.lwd, col=col_lambda[i], lty=2, add=TRUE)
  }
  lambda = obj$lambda0
  curve(lambda, lwd=this.lwd, col="black", lty=1, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Speciation", line=0.75, cex=1.5)


  mu = obj$mu0
  Y_MIN <- min( mu(times) )
  Y_MAX <- max( mu(times) )
  for (i in 1:num.models) {
    mu = obj$mu1[[i]]
    Y_MIN <- min(Y_MIN, mu(times))
    Y_MAX <- max(Y_MAX, mu(times))
  }
  plot( NA, NA, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  for (i in 1:num.models) {
    mu = obj$mu1[[i]]
    curve(mu, lwd=this.lwd, col=col_mu[i], lty=2, add=TRUE)
  }
  mu = obj$mu0
  curve(mu, lwd=this.lwd, col="black", lty=1, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Extinction", line=0.75, cex=1.5)


  lambda = obj$lambda0
  mu     = obj$mu0
  delta  = function(t) lambda(t) - mu(t)
  Y_MIN <- min( delta(times) )
  Y_MAX <- max( delta(times) )
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    mu     = obj$mu1[[i]]
    delta  = function(t) lambda(t) - mu(t)
    Y_MIN <- min(Y_MIN, delta(times))
    Y_MAX <- max(Y_MAX, delta(times))
  }
  plot( NA, NA, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    mu     = obj$mu1[[i]]
    delta  = function(t) lambda(t) - mu(t)
    curve(delta, lwd=this.lwd, col=col_delta[i], lty=2, add=TRUE)
  }
  lambda = obj$lambda0
  mu     = obj$mu0
  delta  = function(t) lambda(t) - mu(t)
  curve(delta, lwd=this.lwd, col="black", lty=1, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Net-diversification", line=0.75, cex=1.5)


  lambda = obj$lambda0
  mu     = obj$mu0
  eps    = function(t) mu(t) / lambda(t)
  Y_MIN <- min( eps(times) )
  Y_MAX <- max( eps(times) )
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    mu     = obj$mu1[[i]]
    eps    = function(t) mu(t) / lambda(t)
    Y_MIN <- min(Y_MIN, eps(times))
    Y_MAX <- max(Y_MAX, eps(times))
  }
  plot( NA, NA, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  for (i in 1:num.models) {
    lambda = obj$lambda1[[i]]
    mu     = obj$mu1[[i]]
    eps    = function(t) mu(t) / lambda(t)
    curve(eps, lwd=this.lwd, col=col_epsilon[i], lty=2, add=TRUE)
  }
  lambda = obj$lambda0
  mu     = obj$mu0
  eps    = function(t) mu(t) / lambda(t)
  curve(eps, lwd=this.lwd, col="black", lty=1, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Relative extinction", line=0.75, cex=1.5)



}
