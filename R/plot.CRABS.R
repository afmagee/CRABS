#' Plots the rate functions including the pulled rates.
#'
#' @param x An object of class "CRABS"
#' @param ... other parameters
#'
#' @return a patchwork object
#' 
#' @export
#' @examples
#' 
#' data(primates_ebd)
#' lambda <- approxfun(primates_ebd$time, primates_ebd$lambda)
#' mu <- approxfun(primates_ebd$time, primates_ebd$mu)
#' times <- seq(0, max(primates_ebd$time), length.out = 500)
#' 
#' model <- create.model(lambda, mu, times = times)
#' 
#' plot(model)
plot.CRABS <- function( x, ... ) {
  df <- model2df(x)
  
  df_lambda <- df %>% filter(grepl("speciation", rate, ignore.case = TRUE))
  df_mu     <- df %>% filter(rate %in% c("Extinction", "Pulled extinction"))
  df_delta  <- df %>% filter(grepl("net-diversification", rate, ignore.case = TRUE))
  df_relext <- df %>% filter(rate == "Relative extinction")
  
  ylim <- range(bind_rows(df_lambda, df_mu)[["value"]])
  
  p1 <- df_lambda %>%
    ggplot(aes(x = Time, y = value, linetype = rate)) +
    theme_classic() +
    geom_line(color = "darkblue") +
    scale_x_reverse() +
    theme(legend.position = c(0.5, 0.6),
          legend.background = element_blank(),
          legend.title = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_linetype_manual(breaks=c("Speciation", "Pulled speciation"), values=c(1,5)) +
    ylab("rate") +
    labs(title = "Speciation") +
    ylim(ylim)
  
  p2 <- df_mu %>%
    ggplot(aes(x = Time, y = value, linetype = rate)) +
    theme_classic() +
    geom_line(color = "orange") +
    scale_x_reverse() +
    theme(legend.position = c(0.5, 0.6),
          legend.background = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          legend.title = element_blank()) +
    scale_linetype_manual(breaks=c("Extinction", "Pulled extinction"), values=c(1,5)) +
    ylab("rate") +
    labs(title = "Extinction") +
    ylim(ylim)
  
  p3 <- df_delta %>%
    ggplot(aes(x = Time, y = value, linetype = rate)) +
    theme_classic() +
    geom_line(color = "purple") +
    scale_x_reverse() +
    theme(legend.position = c(0.5, 0.5),
          plot.title = element_text(hjust = 0.5),
          legend.background = element_blank(),
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

#' Print method for CRABS object
#'
#' @param x and object of class CRABS
#' @param ... other arguments
#'
#' @return nothing
#'
#' @export
#' @examples
#' data(primates_ebd)
#' lambda <- approxfun(primates_ebd$time, primates_ebd$lambda)
#' mu <- approxfun(primates_ebd$time, primates_ebd$mu)
#' times <- seq(0, max(primates_ebd$time), length.out = 500)
#' 
#' model <- create.model(lambda, mu, times = times)
#' 
#' print(model)
print.CRABS <- function(x, ...){
  cat("Piecewise-linear birth-death model\n")
  cat("Knots:", length(x$times), "\n")
  cat("Delta-tau:", x$delta_t, "\n")
  p <- plot.CRABS(x, ...)
  plot(p)
  invisible()
}
