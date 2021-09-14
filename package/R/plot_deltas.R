create_heatmatrix <- function(model, name, threshold){
  times <- model$times
  
  ## \Delta \lambda_{i+1} = \lambda_{i+1} - \lambda_i
  lambdas <- model$lambda(times)
  l_i_plus_one <- tail(model$lambda(times), n = -1)
  l_i <- head(model$lambda(times), n = -1)
  delta_lambda <- (l_i_plus_one - l_i)/l_i_plus_one 
  xtimes <- (head(times, n = -1) + tail(times, n = -1))/2
  
  # is increasing
  inc <- ifelse(delta_lambda > threshold, 1, 0)
  # if decreasing
  dec <- ifelse(delta_lambda < -threshold, -1, 0)
  
  
  z <- factor(inc+dec)
  
  df <- tibble(delta_lambda = delta_lambda,
               direction = z,
               x = xtimes,
               name = name)
  return(df)
}

plotdata <- function(model_set, threshold){
  model_names <- names(model_set)
  l <- lapply(model_names, 
              function(name) create_heatmatrix(model_set[[name]], name, threshold))
  df <- do.call(rbind, l)
  return(df)
}


summary_trends <- function(model_set, threshold = 0.005, return_data = FALSE){
  ##
  df <- plotdata(model_set, threshold)
  
  # Replot rates
  df2_ <- data.frame(lapply(model_set, function(model) model$lambda(model$times)), times = model_set[[1]]$times)
  df2 <- tidyr::gather(df2_, "name", "lambda", -times)
  # Heatmap 
  p1 <- ggplot(df, aes(x, name, fill = direction)) + 
    geom_tile() +
    scale_x_reverse() +
    scale_fill_brewer(palette = "PRGn") +
    theme_bw()
  
  # plot derivative
  p2 <- ggplot(df, aes(x, delta_lambda, color = name)) + 
    geom_line() +
    geom_hline(yintercept = threshold, linetype = "dashed", color = "red") +
    geom_hline(yintercept = -threshold, linetype = "dashed", color = "red") +
    scale_x_reverse() +
    scale_fill_brewer(palette = "PRGn") +
    theme_bw()
  

  p3 <- ggplot(df2, aes(times, lambda, color = name)) +
    geom_line() + 
    scale_x_reverse() +
    theme_bw()
  
  # 
  # grid.newpage()
  # p <- rbind(ggplotGrob(p3), 
  #            ggplotGrob(p2), 
  #            ggplotGrob(p1))
  # grid.draw(p, size = "last")
  # invisible(p)
  p <- cowplot::plot_grid(p3, p2, p1, ncol = 1, greedy = FALSE, align = "vh")
  if(return_data){
    return(list(heatmap_data = df, rate_data = df2))
  }else{
    return(p)  
  }
}