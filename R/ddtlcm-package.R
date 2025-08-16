#' @keywords internal
"_PACKAGE"
#' @seealso
#' \itemize{
#' \item <https://github.com/limengbinggz/ddtlcm> for the source code
#' and system/software requirements to use `ddtlcm` for your data.
#' }

## usethis namespace: start
## usethis namespace: end
NULL

#' @import data.table
#' @import matrixStats
#' @import phylobase
#' @import testthat
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

