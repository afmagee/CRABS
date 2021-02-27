#' Sample custom functions through time.
#'
#' @param num.epochs The number of rates to draw
#' @param lambda0 The rate at present
#' @param v_p_div The pulled diversification rate at all changepioints
#' @param v_ext1 The extinction rate at all changepoints
#' @param delta_t The width of each grid cell
#' @return Extinction rate at all changepoints
#' @export
#' @example 
#' #TODO
sample.rates <- function(num.epochs, lambda0=NULL, rsample=NULL, rsample0=NULL, autocorrelated=FALSE) {
  
  N_SAMPLES = ifelse( is.null(lambda0), num.epochs+1, num.epochs )
  if ( autocorrelated == FALSE ) {
    # we draw a bunch of iid samples
    new_rates = rsample(N_SAMPLES)
    if( is.null(lambda0) == FALSE ) {
      new_rates = c(lambda0, new_rates)
    }
  } else {
    # we draw autocorrelated rates
    if ( is.null(lambda0) == FALSE ) {
      new_rates = c( lambda0 )
    } else if ( is.null(rsample0) == FALSE ) {
      new_rates = c( rsample0() )
    }
    for ( i in 1:num.epochs ) {
      new_rates[i+1] = rsample(new_rates[i])
    }
  }
  
  return (new_rates)
}
