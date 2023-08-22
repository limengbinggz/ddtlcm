#'  Dirichlet diffusion tree-latent class model (DDT-LCM)
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
#' @description
#' `ddtlcm` is designed for
#'  clustering multivariate binary observations over grouped items while (1) leveraging between-cluster similarities guided 
#'  by an unknown tree that is simultaneously estimated, and (2) accounting for varying degrees of 
#'  shrinkage across major item groups. Classes positioned closer on the tree exhibit more similarities a priori. 
#'  The model guards against potential numerical and statistical instability of classical LCMs especially
#'  when classes are weakly separated under small sample sizes. This is achieved by equipping a LCM with 
#'  a DDT process prior on the class profiles, which are the class-conditional response probabilities. The
#'  posterior inference algorithm is based on Metropolis-Hastings algorithm for sampling the tree structure, and
#'  Gibbs sampler with Polya-Gamma augmentation for the LCM parameters.
#' 
#' @name ddtlcm
#' @import data.table
#' @import matrixStats
#' @import phylobase
#' @importFrom BayesLogit rpg
#' @importFrom data.table `:=`
#' @importFrom extraDistr rdirichlet 
#' @importFrom Matrix kronecker Diagonal diag
#' @importFrom methods as
#' @importFrom poLCA poLCA
#' @importFrom stats as.formula cor dist hclust quantile rbinom rgamma rmultinom rnorm runif sd var
#' @importFrom truncnorm rtruncnorm 
# NULL
library(phylobase)
phylobase::phylobase.options(singleton = "ok")

# Make sure data.table knows we know we're using it
.datatable.aware = TRUE

# Prevent R CMD check from complaining about the use of pipe expressions
# standard data.table variables
if (getRversion() >= "2.15.1")
  utils::globalVariables(c(".", ".I", ".N", ".SD"), utils::packageName())

#> NULL






