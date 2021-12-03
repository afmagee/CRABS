#' Plots for congruent birth-death models.
#'
#' Compared to make.congruence.class.plot, this plots individual curves of rates.
#' This allows vizualizing specific patterns, at the cost of either large files or smaller numbers of plotted curves.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param sample.grid Sampled congruent models to plot. Output of sample.congruence.class.
#' @param max.curves Maximum number of curves to plot.
#' @return No value returned, plot is generated.
#' @export
#' @examples
#' #TODO
ACDC.plot.congruence.class.spaghetti <- function(func_spec0, func_ext0, max.t, sample.grid, max.curves=500 ) {
  # recover()

  ## plot settings
  col_lambda    = "darkblue"
  col_mu        = "red3"
  col_delta     = "darkorchid4"
  col_epsilon   = "darkgreen"
  this.lwd      = 2

  ## Here we define some global options
  NUM_TIME_DISCRETIZATIONS = 1000
  NUM_RATE_PLOT_DISCR      = 100
  num.epochs               = length(sample.grid$grid.mu[1,])-1
  num.samples              = length(sample.grid$grid.mu[,1])
  num.rates                = num.epochs + 1

  ## Given the global settings, we can compute some general parameters
  times                    = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
  epoch_times              = (0:num.epochs) / num.epochs * max.t

  if ( is.null( func_spec0 ) == FALSE && is.null( func_ext0 ) == FALSE ) {
    ## construct the net-diversification and relative extinction rate functions
    func_net_div0   <- function(t) {
      func_spec0(t) - func_ext0(t)
    }
    func_rel_ext0   <- function(t) {
      func_ext0(t) / func_spec0(t)
    }

  } else {
    func_net_div0   <- NULL
    func_rel_ext0   <- NULL
  }


  ## extract the rate specific grids from the samples
  grid.mu             = sample.grid$grid.mu
  grid.lambda         = sample.grid$grid.lambda
  grid.net_div        = sample.grid$grid.net_div
  grid.rel_ext        = sample.grid$grid.rel_ext
  grid.delta_mu       = sample.grid$grid.delta_mu
  grid.delta_lambda   = sample.grid$grid.delta_lambda
  grid.delta_net_div  = sample.grid$grid.delta_net_div
  grid.delta_rel_ext  = sample.grid$grid.delta_rel_ext

  ## so we can reset user graphical parameters
  old_par <- par()
  par(mfrow=c(4,2),lend=2,mai=c(0.75,0.75,0.5,0.05),omi=rep(0.01,4))

  ## Cap at 1000 curves per panel, otherwise the resulting plot file size will be too big
  to_plot <- rep(TRUE,num.samples)
  if ( num.samples > max.curves ) {
    to_plot[(max.curves+1):num.samples] <- FALSE
  }

  ## Transparent grey so we can see patterns, could scale to number of plots
  spaghetti_col <- adjustcolor("#666666",alpha.f=0.2)

  ## Make plots!
  Y_MIN <- min( func_spec0(times), grid.lambda )
  Y_MAX <- max( func_spec0(times), grid.lambda )
  matplot(rev(epoch_times),t(grid.lambda[to_plot,]),
          ylim=c(Y_MIN,Y_MAX), lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n", xlab="", ylab="", main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Speciation", line=0.75, cex=1.4)
  if ( !is.null(func_spec0) ) {
    lines(rev(times),func_spec0(times),lwd=this.lwd,col=col_lambda,lty=2)
  }

  matplot(rev(epoch_times),t(grid.delta_lambda[to_plot,]),
          lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n",xlab="",ylab="",main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="delta rate", line=2.25, cex=1.25)
  mtext(side=3, text="Delta speciation", line=0.75, cex=1.4)



  Y_MIN <- min( func_ext0(times), grid.mu )
  Y_MAX <- max( func_ext0(times), grid.mu )
  matplot(rev(epoch_times),t(grid.mu[to_plot,]), ylim=c(Y_MIN,Y_MAX),
          lty=1,lwd=1,type="l", col=spaghetti_col,
          xaxt="n", xlab="", ylab="", main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Extinction", line=0.75, cex=1.4)
  if ( !is.null(func_ext0) ) {
    lines(rev(times),func_ext0(times),lwd=this.lwd,col=col_mu,lty=2)
  }

  matplot(rev(epoch_times),t(grid.delta_mu[to_plot,]),
          lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n",xlab="",ylab="",main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="delta rate", line=2.25, cex=1.25)
  mtext(side=3, text="Delta extinction", line=0.75, cex=1.4)



  Y_MIN <- min( func_net_div0(times), grid.net_div )
  Y_MAX <- max( func_net_div0(times), grid.net_div )
  matplot(rev(epoch_times), t(grid.net_div[to_plot,]), ylim=c(Y_MIN,Y_MAX),
          lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n", xlab="", ylab="", main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Net-diversification", line=0.75, cex=1.4)
  if ( !is.null(func_net_div0) ) {
    lines(rev(times),func_net_div0(times),lwd=this.lwd,col=col_delta,lty=2)
  }

  matplot(rev(epoch_times),t(grid.delta_net_div[to_plot,]),
          lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n",xlab="",ylab="",main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="delta rate", line=2.25, cex=1.25)
  mtext(side=3, text="Delta-net-diversification", line=0.75, cex=1.4)



  Y_MIN <- min( func_rel_ext0(times), grid.rel_ext )
  Y_MAX <- max( func_rel_ext0(times), grid.rel_ext )
  matplot(rev(epoch_times), t(grid.rel_ext[to_plot,]), ylim=c(Y_MIN,Y_MAX),
          lty=1,lwd=1,type="l", col=spaghetti_col,
          xaxt="n", xlab="", ylab="", main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="rate", line=2.25, cex=1.25)
  mtext(side=3, text="Relative extinction", line=0.75, cex=1.4)
  if ( !is.null(func_rel_ext0) ) {
    lines(rev(times),func_rel_ext0(times),lwd=this.lwd,col=col_epsilon,lty=2)
  }

  matplot(rev(epoch_times),t(grid.delta_rel_ext[to_plot,]),
          lty=1,lwd=1,type="l",col=spaghetti_col,
          xaxt="n",xlab="",ylab="",main="")
  axis(1,at=pretty(epoch_times),labels=rev(pretty(epoch_times)))
  mtext(side=1, text="time before present", line=2.5, cex=1.25)
  mtext(side=2, text="delta rate", line=2.25, cex=1.25)
  mtext(side=3, text="Delta-relative-extinction", line=0.75, cex=1.4)

  ## reset user graphical parameters
  par(mfrow=old_par$mfrow,lend=old_par$lend,mai=old_par$mai,omi=old_par$omi)

}
