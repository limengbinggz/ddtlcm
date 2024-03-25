#' Parameters for the HCHS dietary recall data example 
#'
#' A list of five variables containing the food items and parameter estimates obtained from
#'   MCMC chains. The real multivariate binary dataset is not included for privacy.
#'   The variables are as follows:
#'
#' \itemize{
#'   \item item_membership_list. A list of G = 7 elements, where the g-th element contains the 
#'  indices of items belonging to major food group g.
#'   \item item_name_list. A list of G = 7 elements, where the g-th element contains the item
#'    labels of items in major food group g, and the name of the g-th element is the major food 
#'    group label.
#'   \item tree_phylo. The maximum a posterior tree estimate with K = 6 leaves obtained from the 
#'    real HCHS data. Class "phylo".
#'   \item class_probability. A K-vector with entries between 0 and 1. The posterior mean estimate
#'    for class probabilities obtained from the real HCHS data.
#'   \item Sigma_by_group. A G-vector greater than 0. The posterior mean estimate for group-specific 
#'   diffusion variances obtained from the real HCHS data.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name parameter_diet
#' @usage data(parameter_diet)
#' @references \url{https://arxiv.org/abs/2306.04700}
"parameter_diet"


#' Result of fitting DDT-LCM to a semi-synthetic data example
#'
#' This is a "ddtlcm" object obtained from running `ddtlcm_fit` to a semi-synthetic dataset with
#'  1000 posterior samples (for the sake of time). See \code{\link{ddtlcm_fit}} for description of 
#'  the object.
#'
#' @docType data
#' @keywords datasets
#' @name result_diet_1000iters
#' @usage data(result_diet_1000iters)
#' @format A list with 8 elements
"result_diet_1000iters"


#' Synthetic data example
#'
#' This list contains one synthetic data with K = 3 latent classes.
#'   The elements are as follows:
#'
#' \itemize{
#'   \item tree_phylo. A "phylo" tree with K = 3 leaves.
#'   \item class_probability. A K-vector with entries between 0 and 1. 
#'   \item item_membership_list. A list of G = 7 elements, where the g-th element contains the 
#'  indices of items belonging to major food group g.
#'   \item item_name_list. A list of G = 7 elements, where the g-th element contains the item
#'    labels of items in major food group g, and the name of the g-th element is the major food 
#'    group label.
#'   \item Sigma_by_group. A G-vector greater than 0. The group-specific diffusion variances.
#'   \item response_matrix. A binary matrix with N = 100 rows and J = 80 columns. Each row contains
#'    a multivariate binary response vector of a synthetic individual.
#'   \item response_prob. A K by J probability matrix. The k-th row contains the item response 
#'    probabilities of class k.
#'   \item tree_with_parameter. A `phylobase::phylo4d` object. Basically the tree_phylo embedded
#'    with additional logit(response_prob) at the leaf nodes.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name data_synthetic
#' @usage data(data_synthetic)
#' @format A list with 8 elements
"data_synthetic"