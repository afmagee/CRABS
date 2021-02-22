library(shiny)
library(gridExtra)
#library(plotly)
library(ggplot2)
#library(ggplotify)
#library(graphics)

ui <- fluidPage(

  # App title ----
  titlePanel("Congruence class explorer"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(


    # Sidebar panel for inputs ----
    sidebarPanel(

      textAreaInput("lambda0", "Estimated speciation rate. Enter either a constant or a f(t), e.g., 2.5+exp(-0.2*t)", rows = 3),
      textAreaInput("mu0", "Estimated extinction rate", rows = 3),
      numericInput("max_t", "Max age", value = 1, min = 0, max = 100),
      numericInput("max_rate", "Max extinction rate", value = 1, min = 0, max = 100),
      numericInput("n_epochs", "Number of epochs (discretizations)", value = 100, min = 1, max = 10000),
      numericInput("n_samples", "Number of samples from congruent models", value = 2500, min = 1E2, max = 1E6)

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot")

    )

  )

)


# Define server logic required to draw a histogram ----
server <- function(input, output) {


  # This expression that generates a curve is wrapped in a call
  # to renderPlot to indicate that:
  #
  # 1. It is "reactive" and therefore should be automatically
  #    re-executed when inputs (input$lambda) change
  # 2. Its output type is a plot
  output$distPlot <- renderPlot({

    col_lambda = "darkblue"
    col_mu     = "red3"
    col_div    = "green"

    this.lwd   = 5

    ## Here we define some global options
    NUM_TIME_DISCRETIZATIONS = 1000
    MAX_T                    = input$max_t
    MAX_RATE                 = input$max_rate
    NUM_EPOCHS               = input$n_epochs
    NUM_SAMPLES              = input$n_samples
    NUM_RATE_DISCR           = 100
    NUM_RATE_PLOT_DISCR      = 100

    ## Given the global settings, we can compute some general parameters
    times               = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * MAX_T
    epoch_times         = (0:NUM_EPOCHS) / NUM_EPOCHS * MAX_T
    delta_t             = MAX_T / NUM_TIME_DISCRETIZATIONS
    rates               = (0:NUM_RATE_DISCR) / NUM_RATE_DISCR * MAX_RATE

    has_lambda0 = FALSE
    has_mu0     = FALSE
    has_mu1     = FALSE

    if ( is.finite( as.numeric(input$lambda0)) ) {

        has_lambda0 <- TRUE

        lambda0     <- as.numeric(input$lambda0)
        func_spec0  <- function(t) {
            rep(lambda0,length(t))
        }

    } else if ( input$lambda0 != "" ) {

        has_lambda0 <- TRUE

        func_spec0  <- function(t) {
            eval(parse(text=input$lambda0))
        }

    }

    if ( is.finite( as.numeric(input$mu0)) ) {

        has_mu0     <- TRUE

        mu0         <- as.numeric(input$mu0)
        func_ext0  <- function(t) {
            rep(mu0,length(t))
        }

    } else if ( input$mu0 != "" ) {

        has_mu0     <- TRUE

        func_ext0   <- function(t) {
            eval(parse(text=input$mu0))
        }

    }

    if ( has_lambda0 && has_mu0 ) {

        ## create vectors of the speciation and extinction rates at the epoch
        v_spec0     <- func_spec0(times)
        v_ext0      <- func_ext0(times)

        func_net_div0   <- function(t) {
            func_spec0(t) - func_ext0(t)
        }
        func_rel_ext0   <- function(t) {
            func_ext0(t) / func_spec0(t)
        }

        lambda_freqs        <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        mu_freqs            <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        net_div_freqs       <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        rel_ext_freqs       <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        delta_lambda_freqs  <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        delta_mu_freqs      <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        delta_net_div_freqs <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))
        delta_rel_ext_freqs <- array(0,dim = c(NUM_RATE_PLOT_DISCR+1,NUM_EPOCHS+1))

        grid.mu             = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.lambda         = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.net_div        = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.rel_ext        = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.delta_mu       = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.delta_lambda   = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.delta_net_div  = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))
        grid.delta_rel_ext  = array(0,dim = c(NUM_SAMPLES,NUM_EPOCHS+1))

        ###  create the pulled diversification rate
        # compute the derivatives
        l            <- v_spec0[-length(times)]
        l_plus_one   <- v_spec0[-1]
        l_derivative <- (l_plus_one - l) / delta_t
        l_derivative <- c(l_derivative[1],l_derivative)

        # finally, add the 1/lambda * lambda dt to the pulled diversification rate
        v_p_div <- v_spec0 - v_ext0 + 1/v_spec0 * l_derivative

        pb                = txtProgressBar(min = 1, max = NUM_SAMPLES, style = 3)
        for (i in 1:NUM_SAMPLES) {

            mu1         <- sample(x=rates, size=NUM_EPOCHS+1,replace=T)
            func_ext1   <- approxfun(epoch_times,mu1)

            ## create vectors of the speciation and extinction rates at the epoch
            v_ext1      <- func_ext1(times)

            ### compute the new lambda
            v_lambda1    <- c()
            v_lambda1[1] <- func_spec0(0)

            for (j in 2:(NUM_TIME_DISCRETIZATIONS+1)) {

                tmp <- 4*v_lambda1[j-1]*delta_t + (v_p_div[j]*delta_t+v_ext1[j]*delta_t-1)^2
                v_lambda1[j] <- (sqrt(tmp) + v_p_div[j]*delta_t+v_ext1[j]*delta_t-1) / (2*delta_t)

            }

            func_spec1   <- approxfun(times,v_lambda1)
            this_lambda  <- func_spec1(epoch_times)
            this_lambda  <- v_lambda1[(0:NUM_EPOCHS)/(NUM_EPOCHS)*NUM_TIME_DISCRETIZATIONS+1]

            this_mu                 = mu1
            this_net_div            = this_lambda - mu1
            this_rel_ext            = mu1 / this_lambda
            grid.mu[i,]             = mu1
            grid.lambda[i,]         = this_lambda
            grid.net_div[i,]        = this_net_div
            grid.rel_ext[i,]        = this_rel_ext
            grid.delta_mu[i,]       = c(this_mu[-1] - this_mu[-(NUM_EPOCHS+1)], this_mu[NUM_EPOCHS+1] - this_mu[NUM_EPOCHS])
            grid.delta_lambda[i,]   = c(this_lambda[-1] - this_lambda[-(NUM_EPOCHS+1)], this_lambda[NUM_EPOCHS+1] - this_lambda[NUM_EPOCHS])
            grid.delta_net_div[i,]  = c(this_net_div[-1] - this_net_div[-(NUM_EPOCHS+1)], this_net_div[NUM_EPOCHS+1] - this_net_div[NUM_EPOCHS])
            grid.delta_rel_ext[i,]  = c(this_rel_ext[-1] - this_rel_ext[-(NUM_EPOCHS+1)], this_rel_ext[NUM_EPOCHS+1] - this_rel_ext[NUM_EPOCHS])

            if ( sum(is.finite(grid.delta_rel_ext[i,])== FALSE) > 0 ) {
               cat("this_rel_ext:\t\t",this_rel_ext,"\n")
            }
            if ( sum(is.finite(grid.delta_net_div[i,])== FALSE) > 0 ) {
               cat("this_net_div:\t\t",this_net_div,"\n")
            }
            setTxtProgressBar(pb, i)
        }
        close(pb)

        max_mu            = max(grid.mu) * 1.05
        max_lambda        = max(grid.lambda) * 1.05
        max_net_div       = ifelse(max(grid.net_div) < 0, max(grid.net_div) * 0.95, max(grid.net_div) * 1.05)
        min_net_div       = ifelse(min(grid.net_div) < 0, min(grid.net_div) * 1.05, min(grid.net_div) * 0.95)
        max_rel_ext       = max(grid.rel_ext) * 1.05
        max_delta_mu      = ifelse(max(grid.delta_mu) < 0, max(grid.delta_mu) * 0.95, max(grid.delta_mu) * 1.05)
        min_delta_mu      = ifelse(min(grid.delta_mu) < 0, min(grid.delta_mu) * 1.05, min(grid.delta_mu) * 0.95)

        max_delta_lambda  = ifelse(max(grid.delta_lambda) < 0, max(grid.delta_lambda) * 0.95, max(grid.delta_lambda) * 1.05)
        min_delta_lambda  = ifelse(min(grid.delta_lambda) < 0, min(grid.delta_lambda) * 1.05, min(grid.delta_lambda) * 0.95)

        max_delta_net_div = ifelse(max(grid.delta_net_div) < 0, max(grid.delta_net_div) * 0.95, max(grid.delta_net_div) * 1.05)
        min_delta_net_div = ifelse(min(grid.delta_net_div) < 0, min(grid.delta_net_div) * 1.05, min(grid.delta_net_div) * 0.95)

        max_delta_rel_ext = ifelse(max(grid.delta_rel_ext) < 0, max(grid.delta_rel_ext) * 0.95, max(grid.delta_rel_ext) * 1.05)
        min_delta_rel_ext = ifelse(min(grid.delta_rel_ext) < 0, min(grid.delta_rel_ext) * 1.05, min(grid.delta_rel_ext) * 0.95)

        PLOT_RATES_LAMBDA         = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_lambda
        PLOT_RATES_MU             = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_mu
        PLOT_RATES_NET_DIV        = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_net_div-min_net_div) + min_net_div
        PLOT_RATES_REL_EXT        = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * max_rel_ext

        PLOT_RATES_DELTA_LAMBDA   = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_lambda-min_delta_lambda) + min_delta_lambda
        PLOT_RATES_DELTA_MU       = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_mu-min_delta_mu) + min_delta_mu
        PLOT_RATES_DELTA_NET_DIV  = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_net_div-min_delta_net_div) + min_delta_net_div
        PLOT_RATES_DELTA_REL_EXT  = (0:NUM_RATE_PLOT_DISCR) / NUM_RATE_PLOT_DISCR * (max_delta_rel_ext-min_delta_rel_ext) + min_delta_rel_ext


        pb                = txtProgressBar(min = 1, max = NUM_SAMPLES, style = 3)
        for (i in 1:NUM_SAMPLES) {


            this_mu              = grid.mu[i,]
            this_lambda          = grid.lambda[i,]
            this_delta_mu        = grid.delta_mu[i,]
            this_delta_lambda    = grid.delta_lambda[i,]
            this_delta_net_div   = grid.delta_net_div[i,]
            this_delta_rel_ext   = grid.delta_rel_ext[i,]

            for ( j in 1:(NUM_EPOCHS+1) ) {
                xxx = findInterval( this_lambda[j], PLOT_RATES_LAMBDA )
                lambda_freqs[xxx,j] <- lambda_freqs[xxx,j] + 1

                xxx = findInterval( this_mu[j], PLOT_RATES_MU )
                mu_freqs[xxx,j] <- mu_freqs[xxx,j] + 1

                xxx = findInterval( this_lambda[j]-this_mu[j], PLOT_RATES_NET_DIV )
                net_div_freqs[xxx,j] <- net_div_freqs[xxx,j] + 1

                xxx = findInterval( this_mu[j]/this_lambda[j], PLOT_RATES_REL_EXT )
                rel_ext_freqs[xxx,j] <- rel_ext_freqs[xxx,j] + 1

                xxx = findInterval( this_delta_lambda[j], PLOT_RATES_DELTA_LAMBDA )
                delta_lambda_freqs[xxx,j] <- delta_lambda_freqs[xxx,j] + 1

                xxx = findInterval( this_delta_mu[j], PLOT_RATES_DELTA_MU )
                delta_mu_freqs[xxx,j] <- delta_mu_freqs[xxx,j] + 1

                xxx = findInterval( this_delta_net_div[j], PLOT_RATES_DELTA_NET_DIV )
                delta_net_div_freqs[xxx,j] <- delta_net_div_freqs[xxx,j] + 1

                xxx = findInterval( this_delta_rel_ext[j], PLOT_RATES_DELTA_REL_EXT )
                delta_rel_ext_freqs[xxx,j] <- delta_rel_ext_freqs[xxx,j] + 1
            }

            setTxtProgressBar(pb, i)
        }
        close(pb)

        mid_lambda          = max(lambda_freqs) / 2.0
#        mid_lambda          = median(lambda_freqs)
        mid_mu              = max(mu_freqs) / 2.0
        mid_net_div         = max(net_div_freqs) / 2.0
        mid_rel_ext         = max(rel_ext_freqs) / 2.0
        mid_delta_lambda    = max(delta_lambda_freqs) / 2.0
        mid_delta_mu        = max(delta_mu_freqs) / 2.0
        mid_delta_net_div   = max(delta_net_div_freqs) / 2.0
        mid_delta_rel_ext   = max(delta_rel_ext_freqs) / 2.0

        cat("Median:\t\t",mid_lambda,"\n")

        tmp <- expand.grid(PLOT_RATES_LAMBDA, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (lambda_freqs) ) ) )
        p1 <- ggplot(res, aes(time, rate, z = freq)) +
              geom_raster(aes(fill = freq)) +
              scale_fill_gradient2(midpoint = mid_lambda, low="white", mid="orange", high="red2") +
              theme_classic() +
              ggtitle("Speciation") +
              theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
              stat_function(fun=func_spec0, size=2)

        tmp <- expand.grid(PLOT_RATES_DELTA_LAMBDA, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_lambda_freqs) ) ) )
        p2 <- ggplot(res, aes(time, rate, z = freq)) +
              geom_raster(aes(fill = freq)) +
              scale_fill_gradient2(midpoint = mid_delta_lambda, low="white", mid="orange", high="red2") +
              theme_classic() +
              ggtitle("Delta-speciation") +
              theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
#              stat_function(fun=func_spec0, size=2)

        tmp <- expand.grid(PLOT_RATES_MU, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (mu_freqs) ) ) )
        p3 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_mu, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Extinction") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                     stat_function(fun=func_ext0, size=2)

        tmp <- expand.grid(PLOT_RATES_DELTA_MU, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_mu_freqs) ) ) )
        p4 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_delta_mu, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Delta-extinction") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) #+
#                     stat_function(fun=func_ext0, size=2)

        tmp <- expand.grid(PLOT_RATES_NET_DIV, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (net_div_freqs) ) ) )
        p5 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_net_div, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Net-diversification") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                     stat_function(fun=func_net_div0, size=2)

        tmp <- expand.grid(PLOT_RATES_DELTA_NET_DIV, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_net_div_freqs) ) ) )
        p6 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_delta_net_div, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Delta-net-diversification") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) #+
#                     stat_function(fun=func_net_div0, size=2)


        tmp <- expand.grid(PLOT_RATES_REL_EXT, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (rel_ext_freqs) ) ) )
        p7 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_rel_ext, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Relative extinction") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                     stat_function(fun=func_rel_ext0, size=2)

        tmp <- expand.grid(PLOT_RATES_DELTA_REL_EXT, epoch_times)
        res <- as.data.frame(list( rate=tmp[,1], time=tmp[,2], freq=as.numeric( (delta_rel_ext_freqs) ) ) )
        p8 <- ggplot(res, aes(time, rate, z = freq)) +
                     geom_raster(aes(fill = freq)) +
                     scale_fill_gradient2(midpoint = mid_delta_rel_ext, low="white", mid="orange", high="red2") +
                     theme_classic() +
                     ggtitle("Delta-relative extinction") +
                     theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) #+
#                     stat_function(fun=func_rel_ext0, size=2)

        p <- grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, ncol = 2)
#        ggplotly(p)
#        ply1 <- ggplotly(p1)
#        ply2 <- ggplotly(p2)
#        ply3 <- ggplotly(p3)
#        ply4 <- ggplotly(p4)
#        ply5 <- ggplotly(p5)
#        ply6 <- ggplotly(p6)
#        ply7 <- ggplotly(p7)
#        ply8 <- ggplotly(p8)

#        subplot(ply1, ply2, ply3, ply4, ply5, ply6, ply7, ply8, ncol=2)


    }

  }, res=92, height=640)

}



shinyApp(ui = ui, server = server)
