# compute.pulled.diversification <- function( v_spec0, v_ext0, delta_t ) {
#   
#   # compute the derivatives
#   l            <- head(v_spec0, n = -1) #v_spec0[-length(v_ext0)]
#   l_plus_one   <- tail(v_spec0, n = -1) #v_spec0[-1]
#   l_derivative <- (l_plus_one - l) / delta_t
#   l_derivative <- c(l_derivative[1], l_derivative)
#   
#   # finally, add the 1/lambda * lambda dt to the pulled diversification rate
#   v_p_div <- v_spec0 - v_ext0 + (1/v_spec0) * l_derivative
#   
#   return (v_p_div)
# }