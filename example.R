library(ggplot2)
library(gridExtra)
library(TreePar)

# I wasn't feeling authoritative enough to name it myself.
library(PACKAGENAME)

# source("src/utils.R")
# source("src/sample.congruence.class.R")
# source("src/plot.congruence.class.R")
#source("src/nonparametric.pulled.diversification.R")

DATASET = "primates"
tree <- read.tree(file=paste0("data/",DATASET,".tre"))

times <- sort( as.numeric( branching.times( tree ) ) )

x     <- sort(getx(tree),decreasing=TRUE)
max_t <- max(x)


## RevBayes HSMRF example
samples <- read.table(file="output/HSMRFBDP_primates.log",stringsAsFactors=FALSE,header=TRUE)
par <- samples[,grepl(paste0("^","speciation"),names(samples))]
par <- par[,grepl("[0-9]",names(par))]
est_speciation <- apply(par,2,quantile,prob=0.5)
lambda <- approxfun( (0:(length(est_speciation)-1))/(length(est_speciation)-1) * max_t, est_speciation )
par <- samples[,grepl(paste0("^","extinction"),names(samples))]
par <- par[,grepl("[0-9]",names(par))]
est_extinction <- apply(par,2,quantile,prob=0.5)
mu <- approxfun( (0:(length(est_extinction)-1))/(length(est_extinction)-1) * max_t, est_extinction )

ylim <- c(min(0,est_speciation,est_extinction),max(est_speciation,est_extinction))

pdf( paste0(DATASET,"_HSMRF_Spec_Ext.pdf") )
    curve( lambda, xlim=c(0, max_t), lwd=2, lty=1, col="blue", ylim=ylim, main=DATASET, xlab="time", ylab="rate", cex.main=1.5, cex.lab=1.5 )
    curve( mu, add=TRUE, lwd=2, lty=1, col="red" )
    legend("topright",legend=c("speciation","extinction"),col=c("blue","red"), lty=c(1,1), cex=1.2)
dev.off()

est_lambda0 = est_speciation[1]
est_mu0     = est_extinction[1]

rsample_speciation = function(n) runif(n,est_lambda0*1.5,est_lambda0*1.2)
speciation_rate_samples <- function() { sample.rates( num.epochs=100, lambda0=est_lambda0, rsample=rsample_speciation, rsample0=NULL, autocorrelated=FALSE) }

rsample_speciation = function( prev ) rlnorm(n=1,log(prev),0.58/100)
speciation_rate_samples <- function() { sample.rates( num.epochs=100, lambda0=est_lambda0, rsample=rsample_speciation, rsample0=NULL, autocorrelated=TRUE) }

rsample_extinction = function(n) rlnorm(n,log(est_mu0),2*0.58)
extinction_rate_samples <- function() { sample.rates( num.epochs=100, rsample=rsample_extinction, rsample0=NULL, autocorrelated=FALSE) }

rsample_extinction = function( prev ) rlnorm(n=1,log(prev-0.01),0.58/10)
rsample_extinction0 = function() rlnorm(n=1,log(est_mu0),1*0.58)+2

extinction_rate_samples <- function() { sample.rates( num.epochs=100, rsample=rsample_extinction, rsample0=rsample_extinction0, autocorrelated=TRUE) }

samples <- sample.congruence.class(func_spec0=lambda, func_ext0=mu, max.t=max_t, num.epochs=100, num.samples=1000, rate.type="extinction", sample.speciation.rates=speciation_rate_samples, sample.extinction.rates=extinction_rate_samples)
p      = plot.congruence.class(func_spec0=lambda, func_ext0=mu, max.t=max_t, sample.grid=samples )
ggsave(p,file=paste0(DATASET,"_HSMRF_ACDC.pdf"))




## TreePar example
RUN_TREEPAR = FALSE
if ( RUN_TREEPAR ) {
res   <-bd.shifts.optim(x,sampling=c(1,1,1,1),grid=max_t/10,start=0.05*max_t,end=max_t,survival=1,posdiv=TRUE)[[2]]

### test if i shifts explain the tree significantly better than i-1 shifts, here i=1:
best <- 0
for (i in 1:(length(res)-1)) {
    test<-pchisq(2*(res[[i]][1]-res[[i+1]][1]),3)
    #if test>0.95 then i shifts is significantly better than i-1 shifts at a 5% error
    if ( test > 0.95 ) {
        best <- i
    } else {
        break
    }
}


best_result <- res[[best+1]]
speciation <- c()
extinction <- c()
for (i in 1:(best+1)) {
    speciation[i] <- best_result[1+i+(best+1)]/(1-best_result[1+i]); extinction[i] <- best_result[1+i+(best+1)]/(1/best_result[1+i]-1);
}

ylim <- c(min(0,speciation,extinction),max(speciation,extinction))
if ( best > 0 ) {
    shift_times <- best_result[(length(best_result)-best+1):length(best_result)]
    tmp_lambda <- function(x) { index <- findInterval(x, shift_times); return ( speciation[index+1] ) }
    tmp_mu     <- function(x) { index <- findInterval(x, shift_times); return ( extinction[index+1] ) }
    lambda <- Vectorize(tmp_lambda)
    mu     <- Vectorize(tmp_mu)
} else {
    lambda <- function(x) { return ( rep(speciation[1],length(x)) ) }
    mu     <- function(x) { return ( rep(extinction[1],length(x)) ) }
}

pdf( paste0(DATASET,"_TreePar_Spec_Ext.pdf") )
    curve( lambda, xlim=c(0, max_t), lwd=2, lty=1, col="blue", ylim=ylim, main=DATASET, xlab="time", ylab="rate", cex.main=1.5, cex.lab=1.5 )
    curve( mu, add=TRUE, lwd=2, lty=1, col="red" )
    legend("topright",legend=c("speciation","extinction"),col=c("blue","red"), lty=c(1,1), cex=1.2)
dev.off()
samples <- sample.congruence.class(func_spec0=lambda, func_ext0=mu, max.t=max_t, max.rate=2, num.epochs=1000, num.samples=1000)
p      = plot.congruence.class(func_spec0=lambda, func_ext0=mu, max.t=max_t, sample.grid=samples )
ggsave(p,file=paste0(DATASET,"_TreePar_ACDC.pdf"))
}




## non-parametric example
#func_r_p <- nonparametric.pulled.diversification(times)
#
#samples <- sample.congruence.class(func_spec0=NULL, func_ext0=NULL, #func_p_div=func_r_p, max.t=5, max.rate=2, num.epochs=100, num.samples=1000)
#p      = plot.congruence.class(func_spec0=NULL, func_ext0=NULL, max.t=5, #sample.grid=samples )
#print(p)
