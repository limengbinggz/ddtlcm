% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglikehoods.R
\name{logllk_location}
\alias{logllk_location}
\title{Compute log likelihood of parameters}
\usage{
logllk_location(
  tree_phylo4d,
  Sigma_by_group,
  item_membership_list,
  dist_mat = NULL,
  tol = 1e-07
)
}
\arguments{
\item{tree_phylo4d}{a "phylo4d" object}

\item{Sigma_by_group}{a vector of diffusion variances of G groups}

\item{item_membership_list}{a list of G elements, where the g-th element contains the column
indices of \code{data} corresponding to items in major group g}

\item{dist_mat}{a tree-structured covariance matrix from a given tree. Default is NULL. If given
a matrix, then computation of the covariance matrix will be skipped to save time. This is useful
in the Metropolis-Hasting algorithm when the previous proposal is not accepted.}

\item{tol}{a small number to prevent underflow when computing eigenvalues}
}
\value{
A list of two elements: a numeric loglikelihood, a covariance matrix of the input tree
}
\description{
Compute the marginal log likelihood of the parameters on the leaves of a tree
}
\seealso{
Other likelihood functions: 
\code{\link{logllk_ddt}()},
\code{\link{logllk_ddt_lcm}()},
\code{\link{logllk_div_time_one}()},
\code{\link{logllk_div_time_two}()},
\code{\link{logllk_lcm}()},
\code{\link{logllk_tree_topology}()}
}
\concept{likelihood functions}
