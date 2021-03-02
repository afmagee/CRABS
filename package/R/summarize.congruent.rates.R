#' Computes summary statistics for sampled congruent rates.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param sample.grid Sampled congruent models to plot. Output of sample.congruence.class.
#' @return List of summaries for all congruent rates and observed rates.
#' @export
#' @example 
#' #TODO
summarize.congruent.rates <- function(func_spec0,func_ext0,max.t,sample.grid) {
  # recover()
  
  ## Here we define some global options
  NUM_TIME_DISCRETIZATIONS = 1000
  NUM_RATE_PLOT_DISCR      = 100
  num.epochs               = length(sample.grid$grid.mu[1,])-1
  num.samples              = length(sample.grid$grid.mu[,1])
  
  ## Given the global settings, we can compute some general parameters
  times                    = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
  epoch_times              = (0:num.epochs) / num.epochs * max.t
  
  ## Given the global settings, we can compute some general parameters
  times               = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
  epoch_times         = (0:num.epochs) / num.epochs * max.t
  delta_t             = max.t / NUM_TIME_DISCRETIZATIONS
  
  
  ## create vectors of the speciation and extinction rates at the epoch
  v_spec0     <- func_spec0(times)
  v_ext0      <- func_ext0(times)
  v_nd0       <- v_spec0 - v_ext0
  v_re0       <- v_ext0/v_spec0
  
  summary.observed <- rbind(
    c(min(v_spec0),max(v_spec0),span(v_spec0),MASV(v_spec0),FC(v_spec0),fnmean(v_spec0)),
    c(min(v_ext0),max(v_ext0),span(v_ext0),MASV(v_ext0),FC(v_ext0),fnmean(v_ext0)),
    c(min(v_nd0),max(v_nd0),span(v_nd0),MASV(v_nd0),FC(v_nd0),fnmean(v_nd0)),
    c(min(v_re0),max(v_re0),span(v_re0),MASV(v_re0),FC(v_re0),fnmean(v_re0))
  )
  colnames(summary.observed)   <- c("min","max","span","MASV","FC","mean")
  row.names(summary.observed)  <- c("lambda","mu","net_div","relative_ext")

  summary.mu <- cbind(
    apply(sample.grid$grid.mu,1,min),
    apply(sample.grid$grid.mu,1,max),
    apply(sample.grid$grid.mu,1,span),
    apply(sample.grid$grid.mu,1,MASV),
    apply(sample.grid$grid.mu,1,FC),
    apply(sample.grid$grid.mu,1,fnmean)
  )
  colnames(summary.mu) <- c("min","max","span","MASV","FC","mean")

  summary.lambda <- cbind(
    apply(sample.grid$grid.lambda,1,min),
    apply(sample.grid$grid.lambda,1,max),
    apply(sample.grid$grid.lambda,1,span),
    apply(sample.grid$grid.lambda,1,MASV),
    apply(sample.grid$grid.lambda,1,FC),
    apply(sample.grid$grid.lambda,1,fnmean)
  )
  colnames(summary.lambda) <- c("min","max","span","MASV","FC","mean")

  summary.net_div <- cbind(
    apply(sample.grid$grid.net_div,1,min),
    apply(sample.grid$grid.net_div,1,max),
    apply(sample.grid$grid.net_div,1,span),
    apply(sample.grid$grid.net_div,1,MASV),
    apply(sample.grid$grid.net_div,1,FC),
    apply(sample.grid$grid.net_div,1,fnmean)
  )
  colnames(summary.net_div) <- c("min","max","span","MASV","FC","mean")

  summary.rel_ext <- cbind(
    apply(sample.grid$grid.rel_ext,1,min),
    apply(sample.grid$grid.rel_ext,1,max),
    apply(sample.grid$grid.rel_ext,1,span),
    apply(sample.grid$grid.rel_ext,1,MASV),
    apply(sample.grid$grid.rel_ext,1,FC),
    apply(sample.grid$grid.rel_ext,1,fnmean)
  )
  colnames(summary.rel_ext) <- c("min","max","span","MASV","FC","mean")
  
  summaries <- list(summary.mu             = summary.mu,
                    summary.lambda         = summary.lambda,
                    summary.net_div        = summary.net_div,
                    summary.rel_ext        = summary.rel_ext,
                    summary.observed       = summary.observed)
  
  return (summaries)
}

#' Mean absolute sequential variation
#' 
#' Measures net change.
#'
#' @param rates Speciation/extinction rates at beginning/end of each piecewise linear interval.
#' @return MASV
#' @keywords internal
MASV <- function(rates) {
  sum(abs(rates[-1] - rates[-length(rates)]))
}

#' Span of rates
#' 
#' max - min
#'
#' @param rates Speciation/extinction rates at beginning/end of each piecewise linear interval.
#' @return span
#' @keywords internal
span <- function(rates) {
  max(rates) - min(rates)
}

#' Fold change
#' 
#' max(rate(0),rate(t_max)) / min(max(rate(0),rate(t_max)))
#'
#' @param rates Speciation/extinction rates at beginning/end of each piecewise linear interval.
#' @return FC
#' @keywords internal
FC <- function(rates) {
  max(rates[c(1,length(rates))]) / min(rates[c(1,length(rates))])
}

#' Mean value of piecewise linear function
#' 
#' @param rates Speciation/extinction rates at beginning/end of each piecewise linear interval.
#' @return mean
#' @keywords internal
fnmean <- function(rates) {
  w <- 1 # all widths are the same so they might as well be one
  n <- length(rates)
  (sum(rates[1:(n-1)]) + sum(0.5 * (rates[2:n] - rates[1:(n-1)])))/(n-1)
}
