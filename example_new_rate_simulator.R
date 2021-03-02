library(ggplot2)
library(gridExtra)
library(TreePar)

# I wasn't feeling authoritative enough to name it myself.
library(PACKAGENAME)

# source("src/utils.R")
# source("src/sample.congruence.class.R")
# source("src/plot.congruence.class.R")
#source("src/nonparametric.pulled.diversification.R")

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

ylim <- c(min(0,est_speciation,est_extinction),max(est_speciation,est_extinction))

curve( lambda, xlim=c(0, max_t), lwd=2, lty=1, col="blue", ylim=ylim, main="", xlab="time", ylab="rate", cex.main=1.5, cex.lab=1.5 )
curve( mu, add=TRUE, lwd=2, lty=1, col="red" )
legend("topleft",legend=c("speciation","extinction"),col=c("blue","red"), lty=c(1,1), cex=1.2,border=NA,bty="n")

est_lambda0 = est_speciation[1]
est_mu0     = est_extinction[1]

# speciation_rate_samples <- function() { sample.basic.models( num.epochs=100, rate0=est_lambda0, max.rate=2, model="episodic2", fc.mean=1.5) }
# extinction_rate_samples <- function() { sample.basic.models( num.epochs=100, max.rate=2, model="episodic2", fc.mean=1.5) }

speciation_rate_samples <- function() { bootstrap.rates(1000, lambda, max_t, keep.rate0=TRUE, replace=FALSE, block.size=100) }
extinction_rate_samples <- function() { bootstrap.rates(1000, lambda, max_t, keep.rate0=FALSE, reverse.time=TRUE, replace=TRUE, block.size=5) }

samples <- sample.congruence.class(func_spec0=lambda, func_ext0=mu, max.t=max_t, num.epochs=100, num.samples=1000, rate.type="speciation", sample.speciation.rates=speciation_rate_samples, sample.extinction.rates=extinction_rate_samples)
p      = make.congruence.class.plot(func_spec0=lambda, func_ext0=mu, max.t=max_t, sample.grid=samples )
p


