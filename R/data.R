#' Primates phylogenetic tree
#' 
#' The example tree taken from the RevBayes tutorial website
#' 
#' @docType data
#' 
#' @usage data(primates)
#' 
"primates"


#' Primates birth-death model
#' 
#' The results of a bayesian horseshoe markov random field (HSMRF) episodic birth-death model, fitted on the primates tree. One hundred episodes. Each estimate is the posterior median. The time unit is millions of years before the present.
#' 
#' @docType data
#' 
#' @usage data(primates_ebd)
#' 
"primates_ebd"

#' Primates birth-death model
#' 
#' See \code{?primates_ebd}, but including posterior samples instead of a summary.
#' 
#' @docType data
#' 
#' @usage data(primates_ebd_log)
#' 
"primates_ebd_log"