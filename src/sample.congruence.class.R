sample.congruence.class <- function(func_spec0, func_ext0, max.t, max.rate, num.epochs, num.samples) {

    ## Here we define some global options
    NUM_TIME_DISCRETIZATIONS = 1000
    NUM_RATE_DISCR           = 100
    NUM_RATE_PLOT_DISCR      = 100

    ## Given the global settings, we can compute some general parameters
    times               = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
    epoch_times         = (0:num.epochs) / num.epochs * max.t
    delta_t             = max.t / NUM_TIME_DISCRETIZATIONS
    rates               = (0:NUM_RATE_DISCR) / NUM_RATE_DISCR * max.rate


    ## compute the speciation rate at the present
    lambda0     <- func_spec0( 0 )

    ## create vectors of the speciation and extinction rates at the epoch
    v_spec0     <- func_spec0(times)
    v_ext0      <- func_ext0(times)

    grid.mu             = array(0,dim = c(num.samples,num.epochs+1))
    grid.lambda         = array(0,dim = c(num.samples,num.epochs+1))
    grid.net_div        = array(0,dim = c(num.samples,num.epochs+1))
    grid.rel_ext        = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_mu       = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_lambda   = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_net_div  = array(0,dim = c(num.samples,num.epochs+1))
    grid.delta_rel_ext  = array(0,dim = c(num.samples,num.epochs+1))

    ###  create the pulled diversification rate
    v_p_div <- compute.pulled.diversification( v_spec0, v_ext0, times, delta_t )


    pb                = txtProgressBar(min = 1, max = num.samples, style = 3)
    for (i in 1:num.samples) {

        mu1         <- sample(x=rates, size=num.epochs+1,replace=T)
        func_ext1   <- approxfun(epoch_times,mu1)

        ## create vectors of the speciation and extinction rates at the epoch
        v_ext1      <- func_ext1(times)

        ### compute the new lambda
        v_lambda1    <- compute.speciation( lambda0, v_p_div, v_ext1, delta_t )

        func_spec1   <- approxfun(times,v_lambda1)
        this_lambda  <- func_spec1(epoch_times)
        this_lambda  <- v_lambda1[(0:num.epochs)/(num.epochs)*NUM_TIME_DISCRETIZATIONS+1]

        this_mu                 = mu1
        this_net_div            = this_lambda - mu1
        this_rel_ext            = mu1 / this_lambda
        grid.mu[i,]             = mu1
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
