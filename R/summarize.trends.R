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


#' Title
#'
#' @param model_set an object of type "ACDCset"
#' @param threshold a threshold for when \eqn{\Delta\lambdai} should be interpreted as decreasing, flat, or increasing
#' @param rate_name either "lambda" or "mu"
#' @param return_data instead of plots, return the plotting dataframes
#' @param rm_singleton whether or not to remove singeltons. Pass starting at present, going towards ancient
#' @param relative_deltas whether to divide \eqn{\Delta\lambdai} by the local lambda value
#' @param group_names a vector of prefixes, if you want to group the models in a facet. For example 'c("reference", "model")'
#' 
#' @importFrom ggplot2 facet_grid stat
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
#' p <- summarize.trends(model_set, 0.01)
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

  # freq_agree <- df %>% 
  #   group_split(time) %>% 
  #   sapply(agreement)
  
  # delta_times <- df %>% 
  #   dplyr::filter(name == df$name[1]) %>% 
  #   (function(e) e$time)
  
  #df_agree <- tibble::tibble(time = delta_times, freq_agree = freq_agree)
  #
  # p4 <- df_agree %>%
  #   ggplot(aes(x = time, y = freq_agree)) +
  #   geom_col(color = "gray", fill = "gray") +
  #   scale_x_reverse(limits = rev(range(rate_times))) +
  #   theme_bw() +
  #   theme(plot.margin = unit(c(t = 0,r = 0,b = 1,l = 1), "pt"),
  #         axis.title.x=element_blank(),
  #         axis.text.x=element_blank()) +
  #   ylab("Accord") +
  #   ylim(c(0, 1))

  if(return_data){
    return(list(heatmap_data = df
                #rate_data = df2, 
                #df_agree = df_agree
           ))
  }else{
    p <- p1 + p2 + plot_layout(ncol = 1, 
                               #guides = "collect", 
                               heights = c(0.5, 0.5))
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