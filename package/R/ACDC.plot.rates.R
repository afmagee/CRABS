#' Plots the rate functions including the pulled rates.
#'
#' @param obj A list of rate functions produced by .
#' @return void.
#' @export
#' @examples
#' #TODO
ACDC.plot.rates <- function( obj ) {

  ## general settings
  #num.intervals = 1000
  
  #times         = (0:num.intervals) / num.intervals * obj$max.t
  times         <- obj$times#seq(0, to = obj$max.t, by = obj$delta_t)
  num.intervals <- obj$num.intervals

  ## plot settings
  col_lambda    = "darkblue"
  col_mu        = "red3"
  col_delta     = "darkorchid4"
  col_epsilon   = "darkgreen"
  this.lwd      = 3

  ## extract the rate functions
  lambda        = obj$lambda
  #p.lambda      = obj$p.lambda
  # compute the pulled speciation rate
  # this is slow for some reason, so I figured there is no need to compute this before plotting
  #v_p_spec    <- compute.pulled.speciation( sapply(times, lambda), sapply(times, obj$mu), times )
  p.lambda <- pulled.speciation(obj)
  #p.lambda <- approxfun(times,v_p_spec)
  
  mu            = obj$mu
  delta         = obj$delta
  p.delta       = obj$p.delta
  epsilon       = obj$epsilon
  
  
  

  par(mfrow=c(2,2))

  Y_MIN <- min(lambda(times), p.lambda(times))
  Y_MAX <- max(lambda(times), p.lambda(times))
  curve(lambda, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_lambda, ylab="", xlab="", main="", lty=1)
  curve(p.lambda, lwd=this.lwd, col=col_lambda, lty=2, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Speciation", line=0.75, cex=1.4)
  legend("topleft", legend = c("Speciation", "Pulled speciation"), lty = c(1,2), col = col_lambda)

  
  #stop()

  Y_MIN <- min(mu(times))
  Y_MAX <- max(mu(times))
  y <- sapply(times, mu)
  #curve(y, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_mu, ylab="", xlab="", main="", lty=1)
  plot(NULL, xlim=rev(c(0,obj$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  lines(times, y, lwd=this.lwd, col=col_mu, lty=1)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Extinction", line=0.75, cex=1.4)


  Y_MIN <- min(delta(times), p.delta(times))
  Y_MAX <- max(delta(times), p.delta(times))
  #stop()
  curve(delta, xlim=rev(range(times)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_delta, ylab="", xlab="", main="", lty=1)
  curve(p.delta, lwd=this.lwd, col=col_delta, lty=2, add=TRUE)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Net-diversification", line=0.75, cex=1.4)
  legend("topleft", legend = c("Net-diversification", "Pulled net-diversification"), lty = c(1,2), col = col_delta)


  Y_MIN <- min(epsilon(times))
  Y_MAX <- max(epsilon(times))
  curve(epsilon, xlim=rev(range(times)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_epsilon, ylab="", xlab="", main="", lty=1)
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Relative extinction", line=0.75, cex=1.4)

  par(mfrow = c(1,1))

}
