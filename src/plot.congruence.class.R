plot.congruence.class <- function(func_spec0, func_ext0, max.t, sample.grid ) {

    ## Here we define some global options
    NUM_TIME_DISCRETIZATIONS = 1000
    NUM_RATE_PLOT_DISCR      = 100
    num.epochs               = length(sample.grid$grid.mu[1,])-1
    num.samples              = length(sample.grid$grid.mu[,1])

    ## Given the global settings, we can compute some general parameters
    times                    = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
    epoch_times              = (0:num.epochs) / num.epochs * max.t



    if ( is.null( func_spec0 ) == FALSE && is.null( func_ext0 ) == FALSE ) {
        ## construct the net-diversification and relative extinction rate functions
        func_net_div0   <- function(t) {
            func_spec0(t) - func_ext0(t)
        }
        func_rel_ext0   <- function(t) {
            func_ext0(t) / func_spec0(t)
        }

    } else {
        func_net_div0   <- NULL
        func_rel_ext0   <- NULL
    }


    ## extract the rate specific grids from the samples
    grid.mu             = sample.grid$grid.mu
    grid.lambda         = sample.grid$grid.lambda
    grid.net_div        = sample.grid$grid.net_div
    grid.rel_ext        = sample.grid$grid.rel_ext
    grid.delta_mu       = sample.grid$grid.delta_mu
    grid.delta_lambda   = sample.grid$grid.delta_lambda
    grid.delta_net_div  = sample.grid$grid.delta_net_div
    grid.delta_rel_ext  = sample.grid$grid.delta_rel_ext

    lambda_freqs        <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    mu_freqs            <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    net_div_freqs       <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    rel_ext_freqs       <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    delta_lambda_freqs  <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    delta_mu_freqs      <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    delta_net_div_freqs <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))
    delta_rel_ext_freqs <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,num.epochs+1))

    ## compute the max values for each rate type
    max_mu            = max(grid.mu) * 1.05
    max_lambda        = max(grid.lambda) * 1.05
    max_net_div       = ifelse(max(grid.net_div) < 0, max(grid.net_div) * 0.95, max(grid.net_div) * 1.05)
    min_net_div       = ifelse(min(grid.net_div) < 0, min(grid.net_div) * 1.05, min(grid.net_div) * 0.95)
    max_rel_ext       = max(grid.rel_ext) * 1.05
    max_delta_mu      = ifelse(max(grid.delta_mu) < 0, max(grid.delta_mu) * 0.95, max(grid.delta_mu) * 1.05)
    min_delta_mu      = ifelse(min(grid.delta_mu) < 0, min(grid.delta_mu) * 1.05, min(grid.delta_mu) * 0.95)

    ## also compute the max values for the rate deltas
    max_delta_lambda  = ifelse(max(grid.delta_lambda) < 0, max(grid.delta_lambda) * 0.95, max(grid.delta_lambda) * 1.05)
    min_delta_lambda  = ifelse(min(grid.delta_lambda) < 0, min(grid.delta_lambda) * 1.05, min(grid.delta_lambda) * 0.95)

    max_delta_net_div = ifelse(max(grid.delta_net_div) < 0, max(grid.delta_net_div) * 0.95, max(grid.delta_net_div) * 1.05)
    min_delta_net_div = ifelse(min(grid.delta_net_div) < 0, min(grid.delta_net_div) * 1.05, min(grid.delta_net_div) * 0.95)

    max_delta_rel_ext = ifelse(max(grid.delta_rel_ext) < 0, max(grid.delta_rel_ext) * 0.95, max(grid.delta_rel_ext) * 1.05)
    min_delta_rel_ext = ifelse(min(grid.delta_rel_ext) < 0, min(grid.delta_rel_ext) * 1.05, min(grid.delta_rel_ext) * 0.95)

cat("Max mu = ",max_mu,"\n")

    ## now set up a raster of plot values for each rate type
    PLOT_RATES_LAMBDA         = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_lambda
    PLOT_RATES_MU             = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_mu
    PLOT_RATES_NET_DIV        = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_net_div-min_net_div) + min_net_div
    PLOT_RATES_REL_EXT        = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_rel_ext

    PLOT_RATES_DELTA_LAMBDA   = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_lambda-min_delta_lambda) + min_delta_lambda
    PLOT_RATES_DELTA_MU       = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_mu-min_delta_mu) + min_delta_mu
    PLOT_RATES_DELTA_NET_DIV  = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_net_div-min_delta_net_div) + min_delta_net_div
    PLOT_RATES_DELTA_REL_EXT  = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_rel_ext-min_delta_rel_ext) + min_delta_rel_ext


    ## compute the frequencies of observing each values from the samples.
    pb                = txtProgressBar(min = 1, max = num.samples, style = 3)
    for (i in 1:num.samples) {

        this_mu              = grid.mu[i,]
        this_lambda          = grid.lambda[i,]
        this_delta_mu        = grid.delta_mu[i,]
        this_delta_lambda    = grid.delta_lambda[i,]
        this_delta_net_div   = grid.delta_net_div[i,]
        this_delta_rel_ext   = grid.delta_rel_ext[i,]

        for ( j in 1:(num.epochs+1) ) {

            time_index = num.epochs+2-j
            xxx = findInterval( this_lambda[j], PLOT_RATES_LAMBDA )
            lambda_freqs[xxx,time_index] <- lambda_freqs[xxx,time_index] + 1

            xxx = findInterval( this_mu[j], PLOT_RATES_MU )
            mu_freqs[xxx,time_index] <- mu_freqs[xxx,time_index] + 1

            xxx = findInterval( this_lambda[j]-this_mu[j], PLOT_RATES_NET_DIV )
            net_div_freqs[xxx,time_index] <- net_div_freqs[xxx,time_index] + 1

            xxx = findInterval( this_mu[j]/this_lambda[j], PLOT_RATES_REL_EXT )
            rel_ext_freqs[xxx,time_index] <- rel_ext_freqs[xxx,time_index] + 1

            xxx = findInterval( this_delta_lambda[j], PLOT_RATES_DELTA_LAMBDA )
            delta_lambda_freqs[xxx,time_index] <- delta_lambda_freqs[xxx,time_index] + 1

            xxx = findInterval( this_delta_mu[j], PLOT_RATES_DELTA_MU )
            delta_mu_freqs[xxx,time_index] <- delta_mu_freqs[xxx,time_index] + 1

            xxx = findInterval( this_delta_net_div[j], PLOT_RATES_DELTA_NET_DIV )
            delta_net_div_freqs[xxx,time_index] <- delta_net_div_freqs[xxx,time_index] + 1

            xxx = findInterval( this_delta_rel_ext[j], PLOT_RATES_DELTA_REL_EXT )
            delta_rel_ext_freqs[xxx,time_index] <- delta_rel_ext_freqs[xxx,time_index] + 1
        }

        setTxtProgressBar(pb, i)
    }
    close(pb)

    mid_lambda          = max(lambda_freqs) / 2.0
    mid_mu              = max(mu_freqs) / 2.0
    mid_net_div         = max(net_div_freqs) / 2.0
    mid_rel_ext         = max(rel_ext_freqs) / 2.0
    mid_delta_lambda    = max(delta_lambda_freqs) / 2.0
    mid_delta_mu        = max(delta_mu_freqs) / 2.0
    mid_delta_net_div   = max(delta_net_div_freqs) / 2.0
    mid_delta_rel_ext   = max(delta_rel_ext_freqs) / 2.0

    breaks = pretty(epoch_times)

    tmp <- expand.grid(PLOT_RATES_LAMBDA, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (lambda_freqs) ) ) )
    p1 <- ggplot(res, aes(time, rate, z = freq)) +
          geom_raster(aes(fill = freq)) +
          scale_fill_gradient2(midpoint = mid_lambda, low="white", mid="orange", high="red2") +
          theme_classic() +
          ggtitle("Speciation") +
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))
    if ( is.null( func_spec0 ) == FALSE ) {
        p1 <- p1 + stat_function(fun=function(t) func_spec0(max.t-t), size=2)
    }

    tmp <- expand.grid(PLOT_RATES_DELTA_LAMBDA, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_lambda_freqs) ) ) )
    p2 <- ggplot(res, aes(time, rate, z = freq)) +
          geom_raster(aes(fill = freq)) +
          scale_fill_gradient2(midpoint = mid_delta_lambda, low="white", mid="orange", high="red2") +
          theme_classic() +
          ggtitle("Delta-speciation") +
          theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
          scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))

    tmp <- expand.grid(PLOT_RATES_MU, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (mu_freqs) ) ) )
    p3 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_mu, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Extinction") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))
    if ( is.null( func_ext0 ) == FALSE ) {
        p3 <- p3 + stat_function(fun=function(t) func_ext0(max.t-t), size=2)
    }

    tmp <- expand.grid(PLOT_RATES_DELTA_MU, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_mu_freqs) ) ) )
    p4 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_delta_mu, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Delta-extinction") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))

    tmp <- expand.grid(PLOT_RATES_NET_DIV, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (net_div_freqs) ) ) )
    p5 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_net_div, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Net-diversification") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))
    if ( is.null( func_net_div0 ) == FALSE ) {
        p5 <- p5 + stat_function(fun=function(t) func_net_div0(max.t-t), size=2)
    }

    tmp <- expand.grid(PLOT_RATES_DELTA_NET_DIV, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_net_div_freqs) ) ) )
    p6 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_delta_net_div, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Delta-net-diversification") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))


    tmp <- expand.grid(PLOT_RATES_REL_EXT, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (rel_ext_freqs) ) ) )
    p7 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_rel_ext, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Relative extinction") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))
    if ( is.null( func_rel_ext0 ) == FALSE ) {
        p7 <- p7 + stat_function(fun=function(t) func_rel_ext0(max.t-t), size=2)
    }

    tmp <- expand.grid(PLOT_RATES_DELTA_REL_EXT, epoch_times)
    res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_rel_ext_freqs) ) ) )
    p8 <- ggplot(res, aes(time, rate, z = freq)) +
                 geom_raster(aes(fill = freq)) +
                 scale_fill_gradient2(midpoint = mid_delta_rel_ext, low="white", mid="orange", high="red2") +
                 theme_classic() +
                 ggtitle("Delta-relative extinction") +
                 theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                 scale_x_continuous(name ="time before present", breaks=breaks, labels=rev(breaks))

    p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2)

    return (p)
}
