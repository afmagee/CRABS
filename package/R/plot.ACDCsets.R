#' Plots the rate functions
#'
#' @param x A list of objects of class ACDCsets
#' @param ... other parameters
#' @export
#' @examples
#' lambda <- function(x) exp(0.3*x) - 0.5*x + 1
#' mu <- function(x) exp(0.3*x) - 0.2*x + 0.2
#' times <- seq(0, 5, by = 0.005)
#' 
#' model <- create.model(lambda, mu, times = times)
#'
#' mus1 <- list(function(t) 0.2 + exp(0.1*t), 
#'             function(t) 0.2 + sin(0.35*t) + 0.1*t)
#'                         
#' mus2 <- list(function(t) 1.0, 
#'              function(t) 0.5 + 0.2*t)
#'                         
#' models <- list()
#' class(models) <- c("list", "ACDCsets") 
#' models[["set1"]] <- congruent.models(model, mus = mus1)
#' 
#' models[["set2"]] <- congruent.models(model, mus = mus2)
#' 
#' plot(models)
#' 
plot.ACDCsets <- function( x, ... ) {
  op <- graphics::par(mfrow=c(length(x), 4), mar = c(0,3,0,0), oma = c(5,1,3,1))
  rates <- c("lambda", "mu", "delta", "epsilon")
  setnames <- names(x)
  
  # Y_MIN <- Inf
  # Y_MAX <- 0
  # Y_MIN_REL_EXT <- Inf
  # Y_MAX_REL_EXT <- 0
  # 
  # for (models_index in seq_along(x)){
  #   models <- x[[models_index]]
  #   
  #   num.models    = length( models )
  #   
  #   lambda <- models[[1]][["lambda"]]
  #   mu     <- models[[1]][["mu"]]
  #   delta <- models[[1]][["delta"]]
  #   times <- models[[1]][["times"]]
  #   
  #   Y_MIN  <- min( Y_MIN, lambda(times), mu(times), delta(times) )
  #   Y_MAX  <- max( Y_MAX, lambda(times), mu(times), delta(times) )
  # 
  #   for (i in 2:num.models) {
  #     lambda <- models[[i]]$lambda
  #     mu     <- models[[i]]$mu
  #     delta  <- models[[i]]$delta
  #     Y_MIN  <- min(Y_MIN, lambda(times), mu(times), delta(times) )
  #     Y_MAX  <- max(Y_MAX, lambda(times), mu(times), delta(times) )
  #   }
  # 
  #   eps <- models[[1]]$epsilon
  #   Y_MIN_REL_EXT <- min( Y_MIN_REL_EXT, sapply(times, eps) )
  #   Y_MAX_REL_EXT <- max( Y_MAX_REL_EXT, sapply(times, eps) )
  #   for (i in 2:num.models) {
  #     eps <- models[[i]]$epsilon
  #     Y_MIN_REL_EXT <- min(Y_MIN_REL_EXT, sapply(times, eps))
  #     Y_MAX_REL_EXT <- max(Y_MAX_REL_EXT, sapply(times, eps))
  #   }
  # }
  # Y_MINS <- c(Y_MIN, Y_MIN, Y_MIN, Y_MIN_REL_EXT)
  # Y_MAXS <- c(Y_MAX, Y_MAX, Y_MAX, Y_MAX_REL_EXT)
  
  times <- x[[1]][[1]][["times"]]
  Y_MINS <- list()
  Y_MAXS <- list()
  
  for (rate_name in rates){
    Y_MINS[[rate_name]] <- min(sapply(x, function(models) min(sapply(models, function(model) min(model[[rate_name]](times))))))
    Y_MAXS[[rate_name]] <- max(sapply(x, function(models) max(sapply(models, function(model) max(model[[rate_name]](times))))))  
  }
  Y_MINS[["lambda"]] <- Y_MINS[["mu"]] <- do.call(min, Y_MINS[c("lambda", "mu")])
  Y_MAXS[["lambda"]] <- Y_MAXS[["mu"]] <- do.call(max, Y_MAXS[c("lambda", "mu")])

  

  
  for (models_index in seq_along(x)){
    models <- x[[models_index]]
    
    ## general settings
    times <- models[[1]]$times
    num.intervals = length(times)
    max.t <- max(times)
    num.models    = length( models )
    
    ## plot settings
    this.lwd      = 1
    
    # cols = list("lambda" = c("black", tail(brewer.pal(num.models, "Blues"), n = -1)), 
    #             "mu" = c("black", tail(brewer.pal(num.models, "Reds"), n = -1)), 
    #             "delta" = c("black", tail(brewer.pal(num.models, "Purples"), n = -1)), 
    #             "epsilon" = c("black", tail(brewer.pal(num.models, "Greens"), n = -1)))
    cols = list("lambda" = c("black", tail(colorspace::sequential_hcl(palette = "Blues", n = num.models), n = -1)),
                "mu" = c("black", tail(colorspace::sequential_hcl(palette = "Reds", n = num.models), n = -1)),
                "delta" = c("black", tail(colorspace::sequential_hcl(palette = "Purples", n = num.models), n = -1)),
                "epsilon" = c("black", tail(colorspace::sequential_hcl(palette = "Greens", n = num.models), n = -1)))
    

    table1 <- list("lambda" = "Speciation", 
                   "mu" = "Extinction", 
                   "delta" = "Net-diversification", 
                   "epsilon" = "Relative extinction")
    
    for (rate_index in seq_along(rates)) {
      rate_name <- rates[[rate_index]]
      rate <- models[[1]][[rate_name]]

      curve(rate, xlim=rev(c(0,max.t)),
            ylim=c(Y_MINS[[rate_name]],Y_MAXS[[rate_name]]),
            lwd=this.lwd, col="black",
            lty=1, ylab="",
            xlab="", main="", xaxt='n')
      
      for (i in 2:num.models) {
        rate <- models[[i]][[rate_name]]
        lines(times,sapply(times, rate),lwd=this.lwd,col=cols[[rate_name]][i],lty=2)
      }
      
      ## Add annotations, fix axes
      if (models_index == length(x)){
        mtext(side=1, text="time before present", line=2.5, cex=1.25) 
        axis(1)
      }
      else{
        axis(1, labels=FALSE) # add the tick marks but not tick labels
      }
      if (rate_index == 1){
        mtext(side=2, text= setnames[models_index], line=2.25, cex=1.25) 
      }
      if (models_index == 1){
        mtext(side=3, text= table1[[rate_name]], line=0.75, cex=1.5) 
      }
      
    }    
  }
  graphics::par(op)
}



#' Print method for ACDCsets object
#'
#' @param x and object of class ACDCsets
#' @param ... other parameters
#'
#' @export
#' @examples
#' lambda <- function(t) exp(0.3*t) - 0.5*t + 1
#' mu <- function(t) exp(0.3*t) - 0.2*t + 0.2
#' 
#' ## A reference model
#' model <- create.model(lambda, mu, times = seq(0, 5, by = 0.005))
#' 
#' mu1 <- lapply(c(0.5, 1.5, 3.0), function(m) function(t) m)
#' 
#' model_set <- congruent.models(model, mus = mu1)
#' print(model_set)
print.ACDCsets <- function(x, ...){
  cat("Set of piecewise-linear birth-death models\n")
  plot.ACDCsets(x, ...)
  invisible()
}