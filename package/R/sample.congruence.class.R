#' Stochastic exploration of congruent models.
#'
#' @param func_spec0 The speciation rate function (measured in time before present).
#' @param func_ext0 The extinction rate function (measured in time before present).
#' @param func_p_div The pulled diversification rate function (measured in time before present).
#' @param max.t The maximum time (before present) to consider rates.
#' @param num.epochs DESCRIPTION NEEDED
#' @param num.samples DESCRIPTION NEEDED
#' @param rate.type DESCRIPTION NEEDED
#' @param sample.speciation.rates DESCRIPTION NEEDED
#' @param sample.extinction.rates DESCRIPTION NEEDED
#' @return A named list with congruent rates.
#' @export
#' @examples
#' #TODO

sample.congruence.class <- function(func_spec0=NULL, func_ext0=NULL, func_p_div=NULL, max.t, num.epochs, num.samples, rate.type="both", sample.speciation.rates=NULL, sample.extinction.rates=NULL) {

    ## Here we define some global options
    NUM_TIME_DISCRETIZATIONS = 1000
    NUM_RATE_DISCR           = 100
    NUM_RATE_PLOT_DISCR      = 100

    ## Given the global settings, we can compute some general parameters
    times               = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
    epoch_times         = (0:num.epochs) / num.epochs * max.t
    delta_t             = max.t / NUM_TIME_DISCRETIZATIONS


    if ( is.null( func_spec0 ) == FALSE && is.null( func_ext0 ) == FALSE ) {
        ## compute the speciation rate at the present
        lambda0     <- func_spec0( 0 )

        ## create vectors of the speciation and extinction rates at the epoch
        v_spec0     <- func_spec0(times)
        v_ext0      <- func_ext0(times)

        ###  create the pulled diversification rate
        v_p_div <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )
    } else {
        v_p_div <- func_p_div( times )
        lambda0 = 1.0
    }

    grid.mu             = array(0,dim = c(num.samples,num.epochs+1))
    grid.lambda         = array(0,dim = c(num.samples,num.epochs+1))
    grid.net_div        = array(0,dim = c(num.samples,num.epochs+1))
    grid.rel_ext        = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_mu       = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_lambda   = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_net_div  = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_rel_ext  = array(0,dim = c(num.samples,num.epochs+1))



    pb                = txtProgressBar(min = 0, max = num.samples, style = 3)
    setTxtProgressBar(pb, 0)
    for (i in 1:num.samples) {

        if ( rate.type == "extinction" || (rate.type == "both" && (i %% 2) == 1) ) {

            found.valid.sample = FALSE
            while ( found.valid.sample == FALSE ) {
                this_mu     <- sample.rates()

                if ( sum( is.finite( this_mu ) == FALSE ) == 0 ) {
                    found.valid.sample = sum( this_mu < 0 ) == 0
                }
            }
            func_ext1   <- approxfun(epoch_times,this_mu)

            ## create vectors of the speciation and extinction rates at the epoch
            v_ext1      <- func_ext1(times)

            ### compute the new lambda
            v_spec1    <- compute.speciation( lambda0, v_p_div, v_ext1, delta_t )
        } else {

            found.valid.sample = FALSE

            while ( found.valid.sample == FALSE ) {

                this_lambda <- sample.speciation.rates()
                func_spec1  <- approxfun(epoch_times,this_lambda)

                ## create vectors of the speciation and extinction rates at the epoch
                v_spec1     <- func_spec1(times)

                ### compute the new mu
                v_ext1      <- compute.extinction( v_p_div, v_spec1, delta_t )


                if ( sum( is.finite( v_ext1 ) == FALSE ) == 0 ) {
                    found.valid.sample = sum( v_ext1 < 0 ) == 0
                }

#                cat("valid = ",length(found.valid.sample),"\n",sep="")
#                cat("valid = ",ifelse(found.valid.sample,"TRUE","FALSE"),"\n",sep="")
#                cat("mu(",i,") = ",v_ext1,"\n",sep=" ")

            }

#            cat("mu(",i,") = ",v_ext1,"\n",sep="")
        }

        this_lambda            = v_spec1[(0:num.epochs)/(num.epochs)*NUM_TIME_DISCRETIZATIONS+1]
        this_mu                = v_ext1[(0:num.epochs)/(num.epochs)*NUM_TIME_DISCRETIZATIONS+1]

        this_net_div            = this_lambda - this_mu
        this_rel_ext            = this_mu / this_lambda
        grid.mu[i,]             = this_mu
        grid.lambda[i,]         = this_lambda
        grid.net_div[i,]        = this_net_div
        grid.rel_ext[i,]        = this_rel_ext
        grid.delta_mu[i,]       = c(this_mu[-1] - this_mu[-(num.epochs+1)], this_mu[num.epochs+1] - this_mu[num.epochs])
        grid.delta_lambda[i,]   = c(this_lambda[-1] - this_lambda[-(num.epochs+1)], this_lambda[num.epochs+1] - this_lambda[num.epochs])
        grid.delta_net_div[i,]  = c(this_net_div[-1] - this_net_div[-(num.epochs+1)], this_net_div[num.epochs+1] - this_net_div[num.epochs])
        grid.delta_rel_ext[i,]  = c(this_rel_ext[-1] - this_rel_ext[-(num.epochs+1)], this_rel_ext[num.epochs+1] - this_rel_ext[num.epochs])

        setTxtProgressBar(pb, i)
    }
    close(pb)

    res <- list( grid.mu             = grid.mu,
                 grid.lambda         = grid.lambda,
                 grid.net_div        = grid.net_div,
                 grid.rel_ext        = grid.rel_ext,
                 grid.delta_mu       = grid.delta_mu,
                 grid.delta_lambda   = grid.delta_lambda,
                 grid.delta_net_div  = grid.delta_net_div,
                 grid.delta_rel_ext  = grid.delta_rel_ext )

    return (res)
}
