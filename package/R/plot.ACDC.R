#' Plots the rate functions including the pulled rates.
#'
#' @param x An object of class "ACDC"
#' @param ... other parameters
#' 
#' @export
#' @examples
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model <- create.model(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' plot(model)
plot.ACDC <- function( x, ... ) {
  # times         <- x$times
  # num.intervals <- x$num.intervals
  # 
  # ## plot settings
  # col_lambda    = "darkblue"
  # col_mu        = "red3"
  # col_delta     = "darkorchid4"
  # col_epsilon   = "darkgreen"
  # this.lwd      = 3
  # 
  # # compute the pulled speciation rate
  # # this is slow, so I figured there is no need to compute this before plotting
  # p.lambda <- pulled.speciation(x)
  # 
  # ## extract the rate functions
  # lambda        = x$lambda
  # mu            = x$mu
  # delta         = x$delta
  # p.delta       = x$p.delta
  # epsilon       = x$epsilon
  # 
  # op <- par(mfrow=c(2,2), 
  #           mar = c(2,3,3,0), 
  #           oma = c(3,1,0,2)) ## c(bottom, left, top right)
  # 
  # Y_MIN <- min(lambda(times), p.lambda(times))
  # Y_MAX <- max(lambda(times), p.lambda(times))
  # curve(lambda, xlim=rev(c(0,x$max.t)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_lambda, ylab = "", xlab="", main="", lty=1)
  # curve(p.lambda, lwd=this.lwd, col=col_lambda, lty=2, add=TRUE)
  # mtext(side=2, text="rate", line=2.25, cex=1.25)
  # mtext(side=3, text="Speciation", line=0.75, cex=1.4)
  # legend("topleft", legend = c("Speciation", "Pulled speciation"), lty = c(1,2), col = col_lambda)
  # 
  # Y_MIN <- min(mu(times))
  # Y_MAX <- max(mu(times))
  # y <- sapply(times, mu)
  # plot(NULL, xlim=rev(c(0,x$max.t)), ylim=c(Y_MIN,Y_MAX), ylab="", xlab="", main="")
  # lines(times, y, lwd=this.lwd, col=col_mu, lty=1)
  # mtext(side=3, text="Extinction", line=0.75, cex=1.4)
  # 
  # Y_MIN <- min(delta(times), p.delta(times))
  # Y_MAX <- max(delta(times), p.delta(times))
  # curve(delta, xlim=rev(range(times)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_delta, ylab="", xlab="", main="", lty=1)
  # curve(p.delta, lwd=this.lwd, col=col_delta, lty=2, add=TRUE)
  # mtext(side=1, text="time before present", line=2.5, cex=1.25)
  # mtext(side=2, text="rate", line=2.25, cex=1.25)
  # mtext(side=3, text="Net-diversification", line=0.75, cex=1.4)
  # legend("topleft", legend = c("Net-diversification", "Pulled net-diversification"), lty = c(1,2), col = col_delta)
  # 
  # Y_MIN <- min(epsilon(times))
  # Y_MAX <- max(epsilon(times))
  # curve(epsilon, xlim=rev(range(times)), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_epsilon, ylab="", xlab="", main="", lty=1)
  # mtext(side=1, text="time before present", line=2.5, cex=1.25)
  # mtext(side=3, text="Relative extinction", line=0.75, cex=1.4)
  # 
  # par(op) ## Reset the graphical parameters
  
  df <- model2df(x)
  
  df_lambda <- df %>% dplyr::filter(grepl("speciation", rate, ignore.case = TRUE))
  df_mu     <- df %>% dplyr::filter(rate == "Extinction")
  df_delta  <- df %>% dplyr::filter(grepl("net-diversification", rate, ignore.case = TRUE))
  df_relext <- df %>% dplyr::filter(rate == "Relative extinction")
  
  ylim <- range(bind_rows(df_lambda, df_mu)[["value"]])
  
  p1 <- df_lambda %>%
    ggplot(aes(x = Time, y = value, linetype = rate)) +
    theme_classic() +
    geom_line(color = "darkblue") +
    scale_x_reverse() +
    theme(legend.position = c(0.5, 0.6),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_linetype_manual(breaks=c("Speciation", "Pulled speciation"), values=c(1,5)) +
    ylab("rate") +
    labs(title = "Speciation") +
    ylim(ylim)
  
  p2 <- df_mu %>%
    ggplot(aes(x = Time, y = value)) +
    theme_classic() +
    geom_line(color = "red3") +
    scale_x_reverse() +
    theme(legend.position = "NA",
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank()) +
    ylab("rate") +
    labs(title = "Extinction") +
    ylim(ylim)
  
  p3 <- df_delta %>%
    ggplot(aes(x = Time, y = value, linetype = rate)) +
    theme_classic() +
    geom_line(color = "darkorchid4") +
    scale_x_reverse() +
    theme(legend.position = c(0.5, 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank()) +
    scale_linetype_manual(breaks=c("Net-diversification", "Pulled net-diversification"), values=c(1,5)) +
    labs(title = "Net-diversification") +
    ylab("rate") +
    xlab("time before present")
  
  p4 <- df_relext %>%
    ggplot(aes(x = Time, y = value)) +
    theme_classic() +
    geom_line(color = "darkgreen") +
    scale_x_reverse() +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          legend.title = element_blank()) +
    labs(title = "Relative extinction") +
    ylab("rate") +
    xlab("time before present")
  
  p <- p1 + p2 +
    p3 + p4 +
    plot_layout(ncol = 2)
  
  return(p)
}

#' Print method for ACDC object
#'
#' @param x and object of class ACDC
#' @param ... other arguments
#'
#' @export
#' @examples
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' model <- create.model(lambda, mu, times = seq(0, 5, by = 0.005))
#'
#' print(model)
print.ACDC <- function(x, ...){
  cat("Piecewise-linear birth-death model\n")
  cat("Knots:", length(x$times), "\n")
  cat("Delta-tau:", x$delta_t, "\n")
  p <- plot.ACDC(x, ...)
  plot(p)
  invisible()
}
