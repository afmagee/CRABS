#' Summarize trends in the posterior
#'
#' @param posterior a list of CRABS objects, each one representing a sample from the posterior
#' @inheritParams summarize.trends
#' 
#' @return a ggplot object
#' @export summarize.posterior
#' @usage summarize.posterior(posterior, threshold = 0.01, rate_name = "lambda", 
#' return_data = FALSE, rm_singleton = FALSE, per_time = TRUE, 
#' window_size = 1, relative_deltas = FALSE) 
#'
#' @examples
#' data(primates_ebd_log)
#' 
#' posterior <- read.RevBayes(primates_ebd_log, max_t = 65, n_samples = 10)
#' 
#' samples <- sample.congruence.class.posterior(posterior, 
#'                                              num.samples = 5,
#'                                              rate.type = "extinction",
#'                                              rate0.median = 0.1,
#'                                              model = "MRF",
#'                                              max.rate = 1.0)
#' 
#' p <- summarize.posterior(samples, threshold = 0.05)
summarize.posterior <- function(posterior,
                                threshold = 0.01,
                                rate_name = "lambda",
                                return_data = FALSE,
                                rm_singleton = FALSE,
                                per_time = TRUE,
                                window_size = 1,
                                relative_deltas = FALSE){
  times <- posterior[[1]][[1]]$times
  max_t <- max(times)
  n_epochs <- length(times)
  
  res <- list()
  for (i in seq_along(posterior)){
    df <- summarize.trends(posterior[[i]], 
                           threshold = threshold, 
                           window_size = window_size,
                           per_time = per_time,
                           return_data = TRUE)[["heatmap_data"]]
    df["posterior"] <- paste0("sample",i)
    res[[i]] <- df
  }
  plotdata <- bind_rows(res)
  
  if(return_data){
    return(plotdata)
  }
  
  k <- sum(sapply(posterior, length))
  
  #levels_base <- 
  levels1 <- c("-1", "0", "1")
  levels1 <- levels1[levels1 %in% levels(plotdata$direction)]
  
  plotdata$direction <- factor(plotdata$direction, levels = levels1)
  
  p1 <- ggplot(plotdata, aes(x = time, fill = direction)) +
    geom_histogram(aes(y = stat(count / sum(count)*(n_epochs-window_size+1))), binwidth = max_t/n_epochs) +
    theme_classic() +
    scale_fill_manual(values = c("purple","white", "#7fbf7b"), labels = direction_labels) +
    theme(legend.position = c(0.3, 0.3),
          legend.title = element_blank()) +
    scale_x_reverse() +
    labs(y = "model coverage", x = "time before present")
  return(p1)
}