# Filter with increasing regularity thresholds
regularity_filtering <- function(filtering_fraction, samples, samples_penalty){
  eligible_models <- names(samples_penalty[samples_penalty<=quantile(samples_penalty, filtering_fraction)])
  # take 50 samples among the models below the regularity threshold
  n_samples <- min(length(eligible_models), 50)
  filtered_samples <- samples[sample(eligible_models, n_samples)]
  if (!("reference" %in% names(filtered_samples))){
    filtered_samples[[1]] <- samples$reference
    names(filtered_samples)[1] <- "reference"
  }
  class(filtered_samples) <- c("list", "CRABSset")
  return(filtered_samples)
}

plot_rate <- function(samples_list, i, rate, xlab="time before present", 
                      main="Congruent diversification histories"){
  
  palettes = c("lambda" = "Blues", 
               "mu" = "orange", 
               "delta" = "purple", 
               "epsilon" = "Greens")
  
  ylabs = c("lambda" = "Speciation rate", 
            "mu" = "Extinction rate", 
            "delta" = "Net-diversification rate", 
            "epsilon" = "Turnover rate")
  
  sample_i <- samples_list[[i]]
  sample_i <- sample_i[sort(names(sample_i))]
  
  rate_min <- 0.0
  rate_max <- 0.99
  for (i in seq_along(samples_list)){
    samples_set <- samples_list[[i]]
    for (j in seq_along(samples_set)){
      rate_min <- min(rate_min, min(get_rates(samples_set[[j]], rate = rate)))
      rate_max <- max(rate_max, max(get_rates(samples_set[[j]], rate = rate)))
    }
  }
  rate_max <- rate_max * 1.1
  rates_matrix <- sapply(sample_i, get_rates, rate=rate)
  
  cols <- c(head(sequential_hcl(palette=palettes[rate], n=ncol(rates_matrix)),n=-1), "black")
  matplot(sample_i$reference$times, rates_matrix, type="l", lty=c(rep(5,ncol(rates_matrix)-1),1), 
          lwd=c(rep(1,ncol(rates_matrix)-1),1.5), 
          col = cols, ylim=c(rate_min,rate_max), xlim=rev(range(sample_i$reference$times)),
          xlab=xlab, ylab=ylabs[rate], main=main)
}

get_rates <- function(samples, rate) {
  res <- samples[[rate]](samples$times)
  return(res)
}

#' Plots the rate functions after filtering them according to a given penalty and predefined thresholds.
#'
#' @param samples A list of (congruent) CRABS models
#' @param filtering_fractions A vector of thresholds for filtering, as fractions of the most regular trajectories.
#' @param penalty The choice of penalty, among "L1", "L2" and "L1_derivative" (penalty on derivative shifts).
#' @param rates A vector of rate(s) to be plotted, among "lambda" (speciation), "mu" (extinction), "delta" (net-diversification) and "epsilon" (turnover).
#'
#' @return Plots an array of rate trajectories for the chosen rates and thresholds.
#' 
#' @export
#' @examples
#' data("primates_ebd")
#' set.seed(123)
#' 
#' l <- approxfun(primates_ebd[["time"]], primates_ebd[["lambda"]])
#' mu <- approxfun(primates_ebd[["time"]], primates_ebd[["mu"]])
#' times <- primates_ebd[["time"]]
#' 
#' model <- create.model(l, mu, times)
#' 
#' sample.joint.rates <- function(n) {
#'   sample.basic.models.joint(times = times, 
#'                             p.delta = model$p.delta,  
#'                             beta.param = c(0.5,0.3),  
#'                             lambda0 = l(0.0),  
#'                             mu0.median = mu(0.0))
#' }
#' 
#' joint.samples <- sample.congruence.class(model = model, 
#'                                          num.samples = 100, 
#'                                          rate.type = "joint", 
#'                                          sample.joint.rates = sample.joint.rates)
#' 
#' full.plot.regularity.thresholds(joint.samples)
full.plot.regularity.thresholds <- function(samples, 
                                            filtering_fractions=c(0.01, 0.05, 0.2, 0.9), 
                                            penalty="L1", 
                                            rates=c("lambda", "mu")){
  
  # Compute a regularity penalty
  if (penalty=="L1"){
    penalty <- function(model){sum(abs(diff(model$lambda(model$times)))) + sum(abs(diff(model$mu(model$times))))}
  }else if (penalty=="L2"){
    penalty <- function(model){sum(diff(model$lambda(model$times))**2) + sum(diff(model$mu(model$times))**2)}
  }else if (penalty=="L1_derivative"){
    penalty <- function(model){
      return (sum(abs(diff(model$lambda(model$times[-1]))-diff(model$lambda(model$times[-length(model$times)])))) +
                sum(abs(diff(model$mu(model$times[-1]))    -diff(model$mu(model$times[-length(model$times)])))))
    }
  }else{
    stop("Invalid \"penalty\"")
  }
  samples_penalty <- sapply(samples, penalty)

  filtered_samples_list <- lapply(filtering_fractions, regularity_filtering, samples, samples_penalty)

  lr = length(rates)
  lf = length(filtering_fractions)
  par(mfrow=c(lr,lf), mar=c(2.6, 3.6, 2, 0.6), mgp=c(1.5, 0.5, 0))
  for (j in 1:lr){
    for (i in 1:lf){
      main <- ifelse(j==1, paste(filtering_fractions[i]*100, "% most regular trajectories"), "")
      plot_rate(filtered_samples_list, i, rates[j], xlab=ifelse(j==lr, "Time (Mya)", ""), main=main)
    } 
  }
}
