#' Plots the rate functions
#'
#' @param models A list of congruent birth-death models
#' @return nothing
#' @export
#' @examples
#' lambda <- function(x) exp(0.3*x) - 0.5*x + 1
#' mu <- function(x) exp(0.3*x) - 0.2*x + 0.2
#' times <- seq(0, 5, by = 0.005)
#' 
#' model <- create.model(lambda, mu, times = times)
#'
#' mus <- list(function(t) 0.2 + exp(0.1*t), 
#'            function(t) 0.2 + sin(0.35*t) + 0.1*t,
#'                        function(t) 1.0, 
#'                        function(t) 0.5 + 0.2*t)
#' models <- ACDC.congruent.models(model, mus = mus)
#' 
#' plot(models)
plot.ACDCset <- function( models ) {
  
  ## general settings
  times <- models[[1]]$times
  num.intervals = length(times)
  max.t <- max(times)
  num.models    = length( models )
  
  ## plot settings
  col_lambda    = brewer.pal(num.models+1, "Blues")[2:(num.models+1)]
  col_mu        = brewer.pal(num.models+1, "Reds")[2:(num.models+1)]
  col_delta     = brewer.pal(num.models+1, "Purples")[2:(num.models+1)]
  col_epsilon   = brewer.pal(num.models+1, "Greens")[2:(num.models+1)]
  this.lwd      = 1
  
  par(mfrow=c(2,2), mar = c(5, 4, 4, 2) + 0.1, oma = c(0,0,0,0))
  
  lambda <- models[[1]][["lambda"]]
  mu     <- models[[1]][["mu"]]
  Y_MIN  <- min( lambda(times), mu(times) )
  Y_MAX  <- max( lambda(times), mu(times) )
  
  
  
  ## Check 
  for (i in 2:num.models) {
    lambda <- models[[i]]$lambda
    mu     <- models[[i]]$mu
    Y_MIN  <- min(Y_MIN, lambda(times), mu(times) )
    Y_MAX  <- max(Y_MAX, lambda(times), mu(times) )
  }
  lambda = models[[1]]$lambda
  curve(lambda, xlim=rev(c(0,max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col="black", lty=1, ylab="", xlab="", main="")
  for (i in 2:num.models) {
    lambda = models[[i]][["lambda"]]
    lines(times,sapply(times, lambda),lwd=this.lwd,col=col_lambda[i],lty=2)
  }
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Speciation", line=0.75, cex=1.5)
  

  mu <- models[[1]][["mu"]]
  plot(times, sapply(times, mu), "l", xlim=rev(c(0,max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col="black", lty=1, ylab="", xlab="", main="")
  for (i in 2:num.models) {
    mu <- models[[i]][["mu"]]
    lines(times,sapply(times, mu),lwd=this.lwd,col=col_mu[i],lty=2)
  }
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Extinction", line=0.75, cex=1.5)
  
  
  lambda = models[[1]][["lambda"]]
  mu     = models[[1]][["mu"]]
  delta  = function(t) lambda(t) - mu(t)
  Y_MIN <- min( delta(times) )
  Y_MAX <- max( delta(times) )
  for (i in 2:num.models) {
    lambda = models[[i]][["lambda"]]
    mu     = models[[i]][["mu"]]
    delta  = function(t) lambda(t) - mu(t)
    Y_MIN <- min(Y_MIN, delta(times))
    Y_MAX <- max(Y_MAX, delta(times))
  }
  
  lambda = models[[1]][["lambda"]]
  mu     = models[[1]][["mu"]]
  delta  = function(t) lambda(t) - mu(t)
  curve(delta, xlim=rev(c(0,max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col="black", lty=1, ylab="", xlab="", main="")
  for (i in 2:num.models) {
    lambda = models[[i]][["lambda"]]
    mu     = models[[i]][["mu"]]
    delta  = function(t) lambda(t) - mu(t)
    lines(times,delta(times),lwd=this.lwd,col=col_delta[i],lty=2)
  }
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Net-diversification", line=0.75, cex=1.5)
  
  
  lambda = models[[1]][["lambda"]]
  mu     = models[[1]][["mu"]]
  eps    = function(t) mu(t) / lambda(t)
  Y_MIN <- min( sapply(times, eps) )
  Y_MAX <- max( sapply(times, eps) )
  for (i in 2:num.models) {
    lambda = models[[i]][["lambda"]]
    mu     = models[[i]][["mu"]]
    eps    = function(t) mu(t) / lambda(t)
    Y_MIN <- min(Y_MIN, sapply(times, eps))
    Y_MAX <- max(Y_MAX, sapply(times, eps))
  }
  lambda = models[[1]][["lambda"]]
  mu     = models[[1]][["mu"]]
  eps    = function(t) mu(t) / lambda(t)
  curve(eps, xlim=rev(c(0,max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col="black", lty=1, ylab="", xlab="", main="")
  for (i in 2:num.models) {
    lambda = models[[i]][["lambda"]]
    mu     = models[[i]][["mu"]]
    eps    = function(t) mu(t) / lambda(t)
    lines(times,eps(times),lwd=this.lwd,col=col_epsilon[i],lty=2)
  }
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Relative extinction", line=0.75, cex=1.5)
  
  
}