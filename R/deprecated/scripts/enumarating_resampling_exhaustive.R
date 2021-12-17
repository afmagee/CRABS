library(ACDC)
library(TreePar)
library(combinat)

source("R//utils.R")
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
# block permutations
#################

block.size <- 145

rates <- lambda(seq(0,max_t,length.out=NUM_TIME_DISCRETIZATIONS+1))
deltas <- rates[-1] - rates[-(NUM_TIME_DISCRETIZATIONS+1)]

n_blocks <- ceiling(NUM_TIME_DISCRETIZATIONS/block.size)
n_full_blocks <- floor(NUM_TIME_DISCRETIZATIONS/block.size)
blocks <- list()
starts <- (0:(n_full_blocks-1))*block.size
blocks <- lapply(starts,function(i){
  (i+1):(i+block.size)
})
if ( n_blocks != n_full_blocks ) {
  blocks <- c(blocks,list(((n_full_blocks*block.size+1):NUM_TIME_DISCRETIZATIONS)))
}

perms <- permn(1:n_blocks)

#################
# check rate validity
#################

is_valid <- numeric(length(perms))

curves <- matrix(NA,nrow=length(perms),ncol=1001)

for (i in 1:length(perms)) {
  new_blocks <- blocks[perms[[i]]]
  
  resampled_deltas <- unlist(lapply(new_blocks, function(indices){
    deltas[indices]
  }))
  
  this_lambda <- est_lambda0 + c(0,cumsum(resampled_deltas))
  func_spec1  <- approxfun(times,this_lambda)
  
  ## create vectors of the speciation and extinction rates at the epoch
  v_spec1     <- func_spec1(times)
  
  curves[i,] <- v_spec1
  
  ### compute the new mu
  v_ext1      <- compute.extinction( v_p_div, v_spec1, delta_t )
  
  is_valid[i] <- all(v_ext1 > 0)
  
}



