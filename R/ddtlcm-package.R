#'  Dirichlet diffusion tree-latent class model (DDT-LCM)
#'  
#'  `ddtlcm` is designed for
#'  clustering multivariate binary observations over grouped items while (1) leveraging between-cluster similarities guided 
#'  by an unknown tree that is simultaneously estimated, and (2) accounting for varying degrees of 
#'  shrinkage across major item groups. Classes positioned closer on the tree exhibit more similarities a priori. 
#'  The model guards against potential numerical and statistical instability of classical LCMs especially
#'  when classes are weakly separated under small sample sizes. This is achieved by equipping a LCM with 
#'  a DDT process prior on the class profiles, which are the class-conditional response probabilities. The
#'  posterior inference algorithm is based on Metropolis-Hastings algorithm for sampling the tree structure, and
#'  Gibbs sampler with Polya-Gamma augmentation for the LCM parameters.
#'  
#' @seealso
#' \itemize{
#' \item <https://github.com/limengbinggz/ddtlcm> for the source code
#' and system/software requirements to use `ddtlcm` for your data.
#' }
#'
#' @section main ddtlcm wrapper function:
#' [ddtlcm_fit()]
#'
#' @docType package
#' @name ddtlcm
#' @import data.table
#' @import matrixStats
#' @import phylobase
#' @import ape 
#' @importFrom Matrix kronecker Diagonal diag
#' @importFrom extraDistr rdirichlet 
#' @importFrom truncnorm rtruncnorm 
#' @importFrom poLCA poLCA
#' @importFrom BayesLogit rpg
# NULL
library(phylobase)
phylobase::phylobase.options(singleton = "ok")
#> NULL






