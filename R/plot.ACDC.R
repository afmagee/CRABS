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
    geom_line(color = "orange") +
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
    geom_line(color = "purple") +
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
