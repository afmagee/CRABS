#' Plots the rate functions including the pulled rates.
#'
#' @param obj A list of rate functions produced by .
#' @return void.
#' @export
#' @examples
#' #TODO
ACDC.plot.rates <- function( obj ) {

  ## general settings
  num.intervals = 1000
  times         = (0:num.intervals) / num.intervals * obj$max.t

  ## plot settings
  col_lambda    = "darkblue"
  col_mu        = "red3"
  col_delta     = "green"
  col_epsilon   = "purple"
  this.lwd      = 5

  ## extract the rate functions
  lambda        = obj$lambda
  p.lambda      = obj$p.lambda
  mu            = obj$mu
  delta         = obj$delta
  p.delta       = obj$p.delta
  epsilon       = obj$epsilon

  par(mfrow=c(2,2))

  Y_MIN <- min(lambda(times), p.lambda(times))
  Y_MAX <- max(lambda(times), p.lambda(times))
  curve(lambda, xlim=c(0,obj$max.t), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_lambda, ylab="rate", xlab="time before present", lty=1)
  curve(p.lambda, lwd=this.lwd, col=col_lambda, lty=2, add=TRUE)

  Y_MIN <- min(mu(times))
  Y_MAX <- max(mu(times))
  curve(mu, xlim=c(0,obj$max.t), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_mu, ylab="rate", xlab="time before present", lty=1)

  Y_MIN <- min(delta(times), p.delta(times))
  Y_MAX <- max(delta(times), p.delta(times))
  curve(delta, xlim=c(0,obj$max.t), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_delta, ylab="rate", xlab="time before present", lty=1)
  curve(p.delta, lwd=this.lwd, col=col_delta, lty=2, add=TRUE)

  Y_MIN <- min(epsilon(times))
  Y_MAX <- max(epsilon(times))
  curve(epsilon, xlim=c(0,obj$max.t), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_epsilon, ylab="rate", xlab="time before present", lty=1)



}
