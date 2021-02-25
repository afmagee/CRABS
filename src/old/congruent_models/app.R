library(shiny)
library(deSolve)


p.survival <- function(lambda,mu,t) {

    rate <- mu - lambda

    den <- lambda - mu * exp(rate*t)

    res <- (lambda-mu)/ den

}

## This carries out integration for the time-varying
## speciation/extinction model with functions lambda and mu,
## outputting at the times in the vector 'times'.
p.survival.rate <- function(lambda, mu, times) {
  obj <- function(t, state, parms) {
    r <- mu(t) - lambda(t)
    list(c(r))
  }

  tmax <- times[length(times)]
  t.crit <- c(0,tmax)
  n <- 1
  bin <- findInterval(times, t.crit, TRUE)
  y <- c(0)

  xx <- c()
  r_points <- c()
  # we compute the integrals stepwise from t.crit[i] to t.crit[i+1]
  for ( i in seq_len(n) ) {
    j <- i + 1L
    ti <- c(t.crit[[i]], times[bin == i], t.crit[[j]]-1E-9)
    yi <- lsoda(y, ti, obj, tcrit=t.crit[[j]])

    xx <- c(xx,ti)
    r_points <- c(r_points,yi[,2])

    y <- yi[nrow(yi),-1]

  }

  # add a final point slightly over the given time interval. we need this so that some of the numerical routines do not break when calling r(t).
  xx <- c(-1E-4 ,xx, tmax + 1e-8)
  r_points <- c(0,r_points,r_points[length(r_points)])
  rate <- approxfun(xx,r_points)

  # compute the integral int_{x}^{T} mu(t)*exp(r(t)) dt
  f <- function(t) pmin(mu(t)*exp(rate(t)),1E100)
  probs <- array(0,length(xx))
  for (i in length(xx):2) {
    u <- xx[i]
    l <- xx[i-1]
    val <- integrate(f,upper=u,lower=l)$value
    probs[i-1] <- val + probs[i]

  }
  surv <- approxfun(xx,probs)

  list(r=rate,s=surv)

  p_surv <- function( t ) {
    den <- ( 1 + ( ( surv(0) - surv(t) ) / exp(rate(0)) ) )
    return ( 1/den )
  }

  return(p_surv)
}


ui <- fluidPage(

  # App title ----
  titlePanel("Pulled rates and conversion between estimated rate functions"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(


    # Sidebar panel for inputs ----
    sidebarPanel(

      textAreaInput("lambda0", "Estimated speciation rate. Enter either a constant or a f(t), e.g., 2.5+exp(-0.2*t)", rows = 3),
      textAreaInput("mu0", "Estimated extinction rate", rows = 3),
#      textAreaInput("lambda1", "Alternative speciation rate", rows = 3),
      textAreaInput("mu1", "Alternative extinction rate", rows = 3),
      numericInput("max", "Max age", value = 1, min = 0, max = 100),
      numericInput("n", "Number of epochs (discretizations)", value = 100, min = 1, max = 10000)

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

    NUM_EPOCHS = input$n
    MAX_T      = input$max

    times      = (0:NUM_EPOCHS) / NUM_EPOCHS * MAX_T
    delta_t    = MAX_T / NUM_EPOCHS

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

    if ( is.finite( as.numeric(input$mu1)) ) {

        has_mu1     <- TRUE

        mu1         <- as.numeric(input$mu1)
        func_ext1  <- function(t) {
            rep(mu1,length(t))
        }

    } else if ( input$mu1 != "" ) {

        has_mu1     <- TRUE

        func_ext1   <- function(t) {
            eval(parse(text=input$mu1))
        }

    }


    if ( has_lambda0 && has_mu0 ) {

        ## create vectors of the speciation and extinction rates at the epoch
        v_spec0     <- func_spec0(times)
        v_ext0      <- func_ext0(times)


        ###  create the diversification rate function
        func_div0 <- function(t) func_spec0(t) - func_ext0(t)


        ###  create the pulled diversification rate
        # compute the derivatives
        l            <- v_spec0[-length(times)]
        l_plus_one   <- v_spec0[-1]
        l_derivative <- (l_plus_one - l) / delta_t
        l_derivative <- c(l_derivative[1],l_derivative)

        # finally, add the 1/lambda * lambda dt to the pulled diversification rate
        v_p_div <- v_spec0 - v_ext0 + 1/v_spec0 * l_derivative

        func_p_div <- approxfun(times,v_p_div)


        ### create the pulled speciation rate
        #
        p_surv       <- p.survival.rate(func_spec0,func_ext0,times)
#        p_surv       <- function(t) p.survival(func_spec0(0),func_ext0(0),t)
        v_p_spec     <- v_spec0 * p_surv(times)
        func_p_spec <- approxfun(times,v_p_spec)

    }


    if ( has_lambda0 && has_mu0 && has_mu1 ) {

        ## create vectors of the speciation and extinction rates at the epoch
        v_ext1      <- func_ext1(times)

        ### compute the new lambda
        v_lambda1    <- c()
        v_lambda1[1] <- func_spec0(0)

        for (i in 2:(NUM_EPOCHS+1)) {

            tmp <- 4*v_lambda1[i-1]*delta_t + (v_p_div[i]*delta_t+v_ext1[i]*delta_t-1)^2
            v_lambda1[i] <- (sqrt(tmp) + v_p_div[i]*delta_t+v_ext1[i]*delta_t-1) / (2*delta_t)

        }

        func_spec1 <- approxfun(times,v_lambda1)

    }


    Y_MIN <- 0
    Y_MAX <- 0

    if ( has_lambda0 ) {
        Y_MIN <- min(Y_MIN,func_spec0(times))
        Y_MAX <- max(Y_MAX,func_spec0(times))
    }
    if ( has_mu0 ) {
        Y_MIN <- min(Y_MIN,func_ext0(times))
        Y_MAX <- max(Y_MAX,func_ext0(times))
    }
    if ( has_lambda0 && has_mu0 ) {
        Y_MIN <- min(Y_MIN,func_div0(times))
        Y_MIN <- min(Y_MIN,func_p_div(times))
        Y_MIN <- min(Y_MIN,func_p_spec(times))
        Y_MAX <- max(Y_MAX,func_div0(times))
        Y_MAX <- max(Y_MAX,func_p_div(times))
        Y_MAX <- max(Y_MAX,func_p_spec(times))
    }
    if ( has_mu1 ) {
        Y_MIN <- min(Y_MIN,func_ext1(times))
        Y_MAX <- max(Y_MAX,func_ext1(times))
    }
    if ( has_lambda0 && has_mu0 && has_mu1 ) {
        Y_MIN <- min(Y_MIN,func_spec1(times))
        Y_MAX <- max(Y_MAX,func_spec1(times))
    }

    if ( has_lambda0 ) {

        curve(func_spec0, xlim=c(0,MAX_T), ylim=c(Y_MIN,Y_MAX), lwd=this.lwd, col=col_lambda, ylab="rate", xlab="time before present", lty=1)

    }

    if (  has_mu0 ) {

        plot(func_ext0, xlim=c(0,MAX_T), lwd=this.lwd, col=col_mu, lty=1, add=has_lambda0)

    }

    if ( has_lambda0 && has_mu0 ) {

        plot(func_div0, xlim=c(0,MAX_T), lwd=this.lwd, col=col_div, lty=1, add=TRUE)
        points(x=times, y=v_p_spec, lwd=this.lwd, col=col_lambda, pch=4, cex=1.25)
        points(x=times, y=v_p_div, lwd=this.lwd, col=col_div, pch=4, cex=1.25)

    }

    if ( has_mu1 ) {

         plot(func_ext1, xlim=c(0,MAX_T), lwd=this.lwd, col=col_mu, lty=3, add=(has_lambda0||has_mu0))

    }

    if ( has_lambda0 && has_mu0 && has_mu1 ) {

#         plot(x=func_spec1, xlim=c(0,MAX_T), do.points=FALSE, verticals=FALSE, lwd=this.lwd, col=col_lambda, lty=3, add=TRUE)
         points(x=times, y=v_lambda1, lwd=this.lwd, col=col_lambda, pch=20, cex=1.05)
    }

  })

}



shinyApp(ui = ui, server = server)
