% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{create_leaf_cor_matrix}
\alias{create_leaf_cor_matrix}
\title{Create a tree-structured covariance matrix from a given tree}
\usage{
create_leaf_cor_matrix(tree_phylo4d)
}
\arguments{
\item{tree_phylo4d}{a "phylo4d" object}
}
\value{
a K by K covariance matrix
}
\description{
Retrieve the covariance matrix of leaf nodes of a DDT tree
}
\examples{
# load the MAP tree structure obtained from the real HCHS/SOL data
data(data_synthetic)
# extract elements into the global environment
list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
create_leaf_cor_matrix(tree_with_parameter)
}
