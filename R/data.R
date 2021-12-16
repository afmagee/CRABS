#' Primates phylogenetic tree
#' 
#' The example tree taken from the RevBayes tutorial website
#' 
#' @docType data
#' 
#' @usage data(primates)
#' 
"primates"


#' RevBayes Primates birth-death model
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


#' TESS Primates birth-death model
#' 
#' The results of a bayesian episodic birth-death model in the R-package TESS, fitted on the primates tree. One hundred episodes. Each estimate is the posterior median. The time unit is millions of years before the present.
#' 
#' @docType data
#' 
#' @usage data(primates_ebd_tess)
#' 
"primates_ebd_tess"

#' TreePar Primates birth-death model
#' 
#' The results of a birth-death model in the R-package TreePar, fitted on the primates tree. The estimated model has two epochs, that are maximum-likelihood estimates. The time unit is millions of years before the present.
#' 
#' @docType data
#' 
#' @usage data(primates_ebd_treepar)
#' 
"primates_ebd_treepar"