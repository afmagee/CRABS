create_heatmatrix <- function(model, name, rate_name, threshold, relative_deltas = FALSE, gmap = NULL){
  times <- model$times
  delta_t <- model$times[2] - model$times[1]
  
  ## \Delta \lambda_{i+1} = \lambda_{i+1} - \lambda_i
  rate <- model[[rate_name]](times)
  rate_i <- tail(rate, n = -1)
  rate_i_minus_one <- head(rate, n = -1)
  
  if (relative_deltas){
    delta_rate <- (rate_i_minus_one - rate_i) / rate_i
  }else{
    delta_rate <- (rate_i_minus_one - rate_i)
  }
  delta_rate <- delta_rate / delta_t
  
  xtimes <- (head(times, n = -1) + tail(times, n = -1))/2
  
  # is increasing
  inc <- ifelse(delta_rate > threshold, 1, 0)
  # if decreasing
  dec <- ifelse(delta_rate <= -threshold, -1, 0)
  # No change (or flat) is implicitly the value 0
  
  direction <- factor(inc + dec)
  group_name <- gmap[[name]]
  if(is.null(group_name)){
    group_name <- "other"
  }
  
  df <- tibble::tibble(delta_rate = delta_rate,
                       direction = direction,
                       time = xtimes,
                       name = name,
                       group_name = group_name)
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

plotdata <- function(model_set, threshold, rate_name, relative_deltas, gmap){
  model_names <- names(model_set)
  l <- lapply(model_names, 
              function(name) create_heatmatrix(model_set[[name]], name, rate_name, threshold, relative_deltas, gmap))
  df <- do.call(rbind, l)
  
  df$direction <- factor(df$direction, levels = sort(levels(df$direction)))
  return(df)
}


#' Summarize trends in the congruence class
#'
#' @param model_set an object of type "ACDCset"
#' @param threshold a threshold for when \eqn{\Delta \lambda i} should be interpreted as decreasing, flat, or increasing
#' @param rate_name either "lambda" or "mu" or "delta"
#' @param return_data instead of plots, return the plotting dataframes
#' @param rm_singleton whether or not to remove singletons. Pass starting at present, going towards ancient
#' @param relative_deltas whether to divide \eqn{\Delta \lambda i} by the local lambda value
#' @param group_names a vector of prefixes, if you want to group the models in a facet. For example 'c("reference", "model")'
#' 
#' @importFrom ggplot2 facet_grid stat
#'
#' @return a patchwork object
#' @usage summarize.trends(model_set, threshold = 0.005, rate_name = "lambda", 
#' return_data = FALSE, rm_singleton = FALSE, relative_deltas = FALSE, group_names = NULL)
#' @export summarize.trends
#'
#' @examples
#' 
#' data(primates_ebd)
#' lambda <- approxfun(primates_ebd$time, primates_ebd$lambda)
#' mu <- approxfun(primates_ebd$time, primates_ebd$mu)
#' times <- seq(0, max(primates_ebd$time), length.out = 500)
#' 
#' reference <- create.model(lambda, mu, times = times)
#'
#' mus <- list(function(t) exp(0.01*t) - 0.01*t - 0.9,
#'             function(t) exp(-0.02*t) - 0.2,
#'             function(t) exp(-0.07*t) + 0.02*t - 0.5,
#'             function(t) 0.2 + 0.01*t,
#'             function(t) 0.2)
#' 
#' 
#' model_set <- congruent.models(reference, mus = mus)
#' 
#' p <- summarize.trends(model_set, 0.02)
summarize.trends <- function(model_set, 
                           threshold = 0.005, 
                           rate_name = "lambda", 
                           return_data = FALSE, 
                           rm_singleton = FALSE, 
                           relative_deltas = FALSE, 
                           group_names = NULL){
  ##
  if(!is.null(group_names)){
    gmap <- list()
    for (group_name in group_names){
      idx <- which(startsWith(names(model_set), group_name))
      
      if (length(idx) > 0){
        for (i in idx){
          gmap[[names(model_set)[i]]] <- group_name
        }
      }
    }
  }else{
    gmap <- NULL
  }
  
  df <- plotdata(model_set, threshold, rate_name, relative_deltas, gmap)
  
  rate_times <- model_set[[1]]$times
  
  if (rm_singleton){
    l <- df %>% dplyr::group_by(name) %>%
      group_map(~ remove_singletons(.x), .keep = TRUE)
    df <- do.call(rbind, l)
  }
  
  # plot finite-difference derivative
  if(rate_name %in% c("lambda", "mu", "delta")){
    if(relative_deltas){
      lab <- paste0("$\\Delta\\", rate_name, " = \\frac{\\", rate_name, "_{i-1} - \\", rate_name, "_i}{\\Delta t\\", rate_name, "_i}$")
    }else{
      lab <- paste0("$\\Delta\\", rate_name, " = \\frac{\\", rate_name, "_{i-1} - \\", rate_name, "_i}{\\Delta t}$")
    }
    ylabel <- latex2exp::TeX(lab)
  }else{
    stop("rate_name must either be 'lambda' or 'mu' or 'delta'.")
  }
  
  col1 <- list("mu" = "orange",
               "lambda" = "blue",
               "delta" = "purple")[[rate_name]]
  cbPalette <- c(head(sequential_hcl(palette = col1, 
                                     n = length(model_set)), 
                      n = -1), "black")
  
  p1 <- ggplot(df, aes(time, delta_rate, color = name)) + 
    geom_line() +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -threshold, linetype = "dashed", color = "red") +
    scale_x_reverse(limits = rev(range(rate_times))) +
    scale_color_manual(values = cbPalette) +
    theme_classic() +
    ylab(ylabel) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          plot.margin = unit(c(t = 0,r = 1,b = 0,l = 1), "pt"),
          legend.position = "none")

  # plot directions
  p2 <- ggplot(df, aes(time, name, fill = direction)) + 
    geom_tile() +
    scale_x_reverse(limits = rev(range(rate_times))) +
    scale_fill_manual(values = c("purple","white", "#7fbf7b"), labels = direction_labels) +
    theme_classic() +
    theme(plot.margin = unit(c(t = 0,r = 0,b = 1,l = 1), "pt"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank()) +
    ylab("Models") +
    xlab("time before present")
  
  if(!is.null(gmap)){
    p2 <- p2 +
      facet_grid(group_name~., scales="free_y", space="free_y", switch = "y") +
      theme(panel.spacing = unit(c(0), "lines"),
            strip.background = element_rect(fill=NA),
            axis.text.y = element_blank())
      
  }

  if(return_data){
    return(list(heatmap_data = df
           ))
  }else{
    p <- p1 + p2 + plot_layout(ncol = 1, 
                               #guides = "collect", 
                               heights = c(0.5, 0.5))
    return(p)  
  }
}