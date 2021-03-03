library(ACDC)

library(TreePar)

# source("src/utils.R")
# source("src/sample.congruence.class.R")
# source("src/plot.congruence.class.R")
#source("src/nonparametric.pulled.diversification.R")

###################
# empirical rates
###################

tree <- read.nexus("~/Dropbox/UW/Horseshoe_proof_of_concept/2_empirical_analyses/1_Pygopodidae/data/pygo_starting_tree.tre")

times <- sort( as.numeric( branching.times( tree ) ) )

x     <- sort(getx(tree),decreasing=TRUE)
max_t <- max(x)


## RevBayes HSMRF example
samples <- read.table(file="~/Dropbox/UW/Horseshoe_proof_of_concept/2_empirical_analyses/output/HSMRFBDP_10x_run_1.log",stringsAsFactors=FALSE,header=TRUE)
par <- samples[,grepl(paste0("^","speciation"),names(samples))]
par <- par[,grepl("[0-9]",names(par))]
est_speciation <- apply(par,2,quantile,prob=0.5)
lambda <- approxfun( (0:(length(est_speciation)-1))/(length(est_speciation)-1) * max_t, est_speciation )
est_extinction <- rep(quantile(samples$extinction,prob=0.5),length(est_speciation))
mu <- approxfun( (0:(length(est_extinction)-1))/(length(est_extinction)-1) * max_t, est_extinction )

est_lambda0 = est_speciation[1]
est_mu0     = est_extinction[1]

#################
# setup
#################
max.t <- max_t
num.epochs <- 100

## Here we define some global options
NUM_TIME_DISCRETIZATIONS = 1000
NUM_RATE_DISCR           = 100
NUM_RATE_PLOT_DISCR      = 100

## Given the global settings, we can compute some general parameters
times               = (0:NUM_TIME_DISCRETIZATIONS) / NUM_TIME_DISCRETIZATIONS * max.t
epoch_times         = (0:num.epochs) / num.epochs * max.t
delta_t             = max.t / NUM_TIME_DISCRETIZATIONS

###  create the pulled diversification rate
func_spec0 <- lambda
func_ext0 <- mu
v_spec0     <- func_spec0(times)
v_ext0      <- func_ext0(times)

v_p_div <- compute.pulled.diversification( v_spec0, v_ext0, delta_t )


#################
# check rate validity via posterior samples
#################

par <- samples[,grepl(paste0("^","speciation"),names(samples))]
par <- par[,grepl("[0-9]",names(par))]

is_valid <- numeric(dim(samples)[1])

for (i in 1:dim(par)[1]) {
  
  this_speciation <- par[i,]
  func_spec1 <- approxfun( (0:(length(this_speciation)-1))/(length(this_speciation)-1) * max_t, this_speciation )
  
  ## create vectors of the speciation and extinction rates at the epoch
  v_spec1     <- func_spec1(times)
  
  ### compute the new mu
  v_ext1      <- compute.extinction( v_p_div, v_spec1, delta_t )
  
  is_valid[i] <- all(v_ext1 > 0)
  
}
sum(is_valid)


#################
# check rate validity via posterior samples
#################

l1 <- mean(v_spec0[1:190])
l2 <- mean(v_spec0[250:1000])

shift_range <- 150:250

is_valid <- numeric(length(shift_range))

for (i in 1:length(shift_range)) {
  
  shift <- shift_range[i]
  
  tmp <- c(rep(l1,shift),rep(l2,1000-shift+1))
  func_spec1 <- approxfun( (0:1000)/1000 * max.t, tmp )
  
  ## create vectors of the speciation and extinction rates at the epoch
  v_spec1     <- func_spec1(times)
  
  ### compute the new mu
  v_ext1      <- compute.extinction( v_p_div, v_spec1, delta_t )
  
  is_valid[i] <- all(v_ext1 > 0)
  
}



