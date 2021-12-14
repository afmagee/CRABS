#' Plots the rate functions
#'
#' @param x A list of congruent birth-death x
#' @param ... other parameters
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
#' models <- congruent.models(model, mus = mus)
#' 
#' plot(models)
plot.ACDCset <- function( x, ... ) {
  dfs <- lapply(x, model2df)
  ## Add names columns
  for (i in seq_along(dfs)){
    df <- dfs[[i]]
    df$name <- names(dfs)[i]
    dfs[[i]] <- df
  }
  df <- bind_rows(dfs)
  
  df_lambda <- df %>% filter(rate == "Speciation")
  df_mu     <- df %>% filter(rate == "Extinction")
  df_delta  <- df %>% filter(rate == "Net-diversification")
  df_relext <- df %>% filter(rate == "Relative extinction")
  
  ylim <- range(bind_rows(df_lambda, df_mu)[["value"]])
  
  ## Speciation rate
  col_lambda <- c(head(colorspace::sequential_hcl(palette = "Blues", n = length(unique(df_lambda$name))), n = -1), "black")
  p1 <- df_lambda %>%
    ggplot(aes(x = Time, y = value, color = name)) +
    scale_x_reverse() +
    theme_classic() +
    geom_line(data=subset(df_lambda, name == "reference"), linetype=1) +
    geom_line(data=subset(df_lambda, name != "reference"), linetype="longdash") +
    labs(title = "Speciation") +
    theme(legend.position = "NA",
          axis.title.x = element_blank(),
    ) +
    ylim(ylim) +
    ylab("rate") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = col_lambda)
  
  ## Extinction rate
  col_mu <- c(head(colorspace::sequential_hcl(palette = "orange", n = length(unique(df_mu$name))), n = -1), "black")
  p2 <- df_mu %>%
    ggplot(aes(x = Time, y = value, color = name)) +
    scale_x_reverse() +
    theme_classic() +
    geom_line(data=subset(df_mu, name == "reference"), linetype=1) +
    geom_line(data=subset(df_mu, name != "reference"), linetype="longdash") +
    ggtitle("Extinction") +
    ylim(ylim) +
    theme(legend.position = "NA",
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
    ) +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = col_mu)
  
  ## Net-diversification rate
  col_delta <- c(head(colorspace::sequential_hcl(palette = "purple", n = length(unique(df_delta$name))), n = -1), "black")
  p3 <- df_delta %>%
    ggplot(aes(x = Time, y = value, color = name)) +
    scale_x_reverse() +
    theme_classic() +
    geom_line(data=subset(df_delta, name == "reference"), linetype=1) +
    geom_line(data=subset(df_delta, name != "reference"), linetype="longdash") +
    ggtitle("Net-diversification") +
    theme(legend.position = "NA",
    ) +
    ylab("rate") +
    xlab("time before present") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = col_delta)
  
  ## Relative-extinction rate
  col_relext <- c(head(colorspace::sequential_hcl(palette = "green", n = length(unique(df_relext$name))), n = -1), "black")
  p4 <- df_delta %>%
    ggplot(aes(x = Time, y = value, color = name)) +
    scale_x_reverse() +
    theme_classic() +
    geom_line(data=subset(df_relext, name == "reference"), linetype=1) +
    geom_line(data=subset(df_relext, name != "reference"), linetype="longdash") +
    ggtitle("Relative extinction") +
    theme(legend.position = "NA",
          axis.title.y = element_blank(),
    ) +
    xlab("time before present") +
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = col_relext)
  
  p <- p1 + p2 + 
    p3 + p4 +
    plot_layout(ncol = 2)
  
  return(p)
}

#' Print method for ACDCset object
#'
#' @param x an object of class ACDCset
#' @param ... other arguments
#'
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
#' models <- congruent.models(model, mus = mus)
#' 
#' print(models)
print.ACDCset <- function(x, ...){
  cat("A congruent set of piecewise-linear birth-death models\n")
  cat("Knots:", length(x[[1]]$times), "\n")
  cat("Delta-tau:", x[[1]]$delta_t, "\n")
  cat("n_models: ", length(x), "\n")
  if (length(x) <= 50){
    p <- plot.ACDCset(x)  
    plot(p)
  }else{
    cat("Your set is too large (>50), and won't be plotted.")
  }
  
  invisible()
}

