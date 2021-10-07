create_heatmatrix <- function(model, name, rate_name, threshold, relative_deltas = FALSE){
  times <- model$times
  
  ## \Delta \lambda_{i+1} = \lambda_{i+1} - \lambda_i
  rate <- model[[rate_name]](times)
  rate_i_plus_one <- tail(rate, n = -1)
  rate_i <- head(rate, n = -1)
  
  if (relative_deltas){
    delta_rate <- (rate_i - rate_i_plus_one) / rate_i_plus_one
  }else{
    delta_rate <- rate_i - rate_i_plus_one
  }
  
  xtimes <- (head(times, n = -1) + tail(times, n = -1))/2
  
  # is increasing
  inc <- ifelse(delta_rate > threshold, 1, 0)
  # if decreasing
  dec <- ifelse(delta_rate <= -threshold, -1, 0)
  # No change (or flat) is implicitly the value 0
  
  direction <- factor(inc + dec)
  
  df <- tibble::tibble(delta_rate = delta_rate,
                       direction = direction,
                       time = xtimes,
                       name = name)
  return(df)
}

remove_singletons <- function(df){
  for (i in 2:(nrow(df)-1)){
    triplet <- df[["direction"]][(i-1):(i+1)] ## Create the 3-window
    if (triplet[1] == triplet[3]){ ## Check if ends are equal
      if (triplet[2] != triplet[1]){ 
        df[i, "direction"] <- triplet[3] ## If middle is unequal, assign it to the corner values
      }   
    }
  }
  return(df)
}

direction_labels <- function(x){
  ifelse(x == "-1", "Decreasing", ifelse(x == "1", "Increasing", "Flat"))
}

plotdata <- function(model_set, threshold, rate_name, relative_deltas){
  model_names <- names(model_set)
  l <- lapply(model_names, 
              function(name) create_heatmatrix(model_set[[name]], name, rate_name, threshold, relative_deltas))
  df <- do.call(rbind, l)
  return(df)
}


#' Title
#'
#' @param model_set an object of type "ACDCset"
#' @param threshold a threshold for when \eqn{\Delta\lambdai} should be interpreted as decreasing, flat, or increasing
#' @param rate_name either "lambda" or "mu"
#' @param return_data instead of plots, return the plotting dataframes
#' @param rm_singleton whether or not to remove singeltons. Pass starting at present, going towards ancient
#' @param relative_deltas whether to divide \eqn{\Delta\lambdai} by the local lambda value
#'
#' @return a patchwork object
#' @export
#'
#' @examples
#' 
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' reference <- create.model(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' mu1 <- list(function(t) exp(0.2*t) - 0.3*t + 0.4,
#' function(t) exp(-0.8*t) - 0.1*t + 0.8,
#' function(t) exp(-1.5*t) + 0.2*t + 1.2,
#' function(t) 1.2 + 0.5*t,
#' function(t) 1.5 - 0.28*t)
#' 
#' 
#' model_set <- congruent.models(reference, mus = mu1)
#' 
#' p <- summary_trends(model_set, 0.01)
#' 
summary_trends <- function(model_set, threshold = 0.005, rate_name = "lambda", return_data = FALSE, rm_singleton = TRUE, relative_deltas = FALSE){
  ##
  df <- plotdata(model_set, threshold, rate_name, relative_deltas)
  
  rate_times <- model_set[[1]]$times
  
  if (rm_singleton){
    l <- df %>% dplyr::group_by(name) %>%
      group_map(~ remove_singletons(.x), .keep = TRUE)
    df <- do.call(rbind, l)
  }
  
  cbPalette <- c(head(colorspace::sequential_hcl(palette = "blues", n = length(model_set)), n = -1), "black")
  
  # Replot rates
  df2_ <- data.frame(lapply(model_set, function(model) model[[rate_name]](model$times)), times = model_set[[1]]$times)
  #df2 <- tidyr::gather(df2_, "name", "rate_name", -times)
  #browser()
  df2 <- tidyr::pivot_longer(df2_, cols = -times, names_to = "name", values_to = "rate")

  # plot rates
  #p1 <- ggplot(df2, aes(times, lambda, color = name)) +
  p1 <- ggplot(df2, aes_string("times", "rate", color = "name")) +
    geom_line() + 
    scale_x_reverse(limits = rev(range(rate_times))) +
    scale_color_manual(values = cbPalette) +
    theme_bw() +
    ylab(latex2exp::TeX("$\\lambda$")) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          plot.margin = grid::unit(c(t = 1,r = 1,b = 0,l = 1), "pt"))
  
  # plot finite-difference derivative
  if(rate_name == "lambda"){
    if(relative_deltas){
      ylabel <- latex2exp::TeX("$\\Delta\\lambda = \\frac{\\lambda_{i+1} - \\lambda_i}{\\lambda_i}$")
    }else{
      ylabel <- latex2exp::TeX("$\\Delta\\lambda = \\lambda_{i+1} - \\lambda_i$")
    }
  }else if(rate_name == "mu"){
    if(relative_deltas){
      ylabel <- latex2exp::TeX("$\\Delta\\mu = \\frac{\\mu_{i+1} - \\mu}{\\mu}$")
    }else{
      ylabel <- latex2exp::TeX("$\\Delta\\mu = \\mu_{i+1} - \\mu$")
    }
  }else{
    stop("rate_name must either be 'lambda' or 'mu'.")
  }

  
  p2 <- ggplot(df, aes(time, delta_rate, color = name)) + 
    geom_line() +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -threshold, linetype = "dashed", color = "red") +
    scale_x_reverse(limits = rev(range(rate_times))) +
    scale_color_manual(values = cbPalette) +
    theme_bw() +
    ylab(ylabel) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          plot.margin = unit(c(t = 0,r = 1,b = 0,l = 1), "pt"))
  
  # plot directions
  p3 <- ggplot(df, aes(time, name, fill = direction)) + 
    geom_tile() +
    scale_x_reverse(limits = rev(range(rate_times))) +
    scale_fill_manual(values = c("gray", "green","purple"), labels = direction_labels)
    theme_bw() +
    theme(plot.margin = unit(c(t = 0,r = 0,b = 1,l = 1), "pt")) +
    ylab("Direction") +
    xlab("Time before present")

  freq_agree <- df %>% 
    group_split(time) %>% 
    sapply(agreement)
  
  delta_times <- df %>% 
    dplyr::filter(name == "reference") %>% 
    (function(e) e$time)
  
  df_agree <- tibble::tibble(time = delta_times, freq_agree = freq_agree)
  
  p4 <- df_agree %>%
    ggplot(aes(x = time, y = freq_agree)) +
    geom_col(color = "gray", fill = "gray") +
    scale_x_reverse(limits = rev(range(rate_times))) +
    theme_bw() +
    theme(plot.margin = unit(c(t = 0,r = 0,b = 1,l = 1), "pt"),
          axis.title.x=element_blank(),
          axis.text.x=element_blank()) +
    ylab("Accord") +
    ylim(c(0, 1))

  if(return_data){
    return(list(heatmap_data = df, rate_data = df2, df_agree = df_agree))
  }else{
    p <- p1 + p2 + p4 + p3 + plot_layout(ncol = 1, guides = "collect", heights = c(0.3, 0.3, 0.08, 0.32))
    return(p)  
  }
}

agreement <- function(df){
  d <- df[["direction"]]
  #ref <- d[df[["name"]] == "reference"]
  #n <- sum(d == ref) / length(d) ## Agreement with the reference
  n <- max(table(d))/length(d) ## Majority consensus
  return(n)
}