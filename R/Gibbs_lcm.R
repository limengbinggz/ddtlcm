###############################################################
# Gibbs sampler with polya-gamma and logistic data augmentation
# to sample LCM parameters.
###############################################################

#'Sample item group-specific variances through Gibbs sampler
#'@param shape0 a vector of G elements, each being the shape of the
#'    Inverse-Gamma prior of group g
#'@param rate0 a vector of G elements, each being the rate of the
#'    Inverse-Gamma prior of group g
#'@param dist_mat a list, containing the KxK tree-structured matrix of leaf nodes,
#'    where K is the number of leaves / latent classes, and SVD components
#'@param item_membership_list a vector of G elements, each indicating the number of
#'    items in this group
#'@param locations a KxJ matrix of leaf parameters
#'@return a numeric vector of G elements, each being the newly sampled variance
#'    of the latent location of this group
sample_sigmasq <- function(shape0, rate0, dist_mat, item_membership_list, locations){
  # number of item groups
  G <- length(item_membership_list)
  # number of leaves
  K <- nrow(dist_mat$dist_mat_g)
  Sigma_by_group <- rep(0, G)
  for (g in 1:G) {
    shape <- shape0[g] + K * length(item_membership_list[[g]]) * 0.5
    # get the location of the root node
    locations_g_diff <- locations[, item_membership_list[[g]], drop=F]
    # trace of A^t B A = vec(B^t A)^t vec(A)
    rate <- rate0[g] + 0.5 * sum( (dist_mat$Dhalf_Ut %*% locations_g_diff) **2)
    Sigma_by_group[g] <- 1.0 / rgamma(1, shape = shape, rate = rate)
  }
  return(Sigma_by_group)
}



#'Sample the leaf locations and Polya-Gamma auxilliary variables
#'@param item_membership_list a vector of G elements, each indicating the number of
#'    items in this group
#'@param dist_mat_old a list of leaf covariance matrix from the previous iteration. The list
#'    has length G, the number of item groups
#'@param Sigma_by_group a vector of length G, each denoting the variance of the
#'    brownian motion
#'@param pg_mat a K by J matrix of PG variables from the previous iteration
#'@param a_pg a N by J matrix of hyperparameters of the generalized logistic distribution
#'@param auxiliary_mat a N by J matrix of truncated normal variables from previous iteration
#'@param auxiliary_mat_range a list of two named elements: lb and ub. Each is an N by J 
#'  matrix of the lower/upper bounds of the truncated normal variables.
#'@param class_assignments an integer vector of length N for the individual class assignments.
#'  Each element takes value in 1, ..., K.
#'@return a named list of three matrices: the newly sampled leaf parameters, the Polya-gamma random variables,
#'  and the auxiliary truncated normal variables
sample_leaf_locations_pg <- function(item_membership_list, dist_mat_old,
                                     Sigma_by_group, pg_mat, a_pg, auxiliary_mat, 
                                     auxiliary_mat_range, class_assignments){
  N <- nrow(pg_mat)
  # number of classes
  K <- ncol(dist_mat_old$Dhalf_Ut)
  # total number of items
  J <- length(unlist(item_membership_list))
  # cumulative index of item groups
  new_leaf_data <- matrix(NA, nrow = K, ncol = J)
  for (g in seq_along(item_membership_list)) {
    J_g <- length(item_membership_list[[g]])
    precision_mat <- crossprod(dist_mat_old$Dhalf_Ut) / Sigma_by_group[g]
    # extract indices
    indices <- item_membership_list[[g]]
    pg_mat_g <- pg_mat[, indices, drop = FALSE]
    auxiliary_mat_g <- auxiliary_mat[, indices, drop = FALSE]
    xi_1_raw <- pg_mat_g * auxiliary_mat_g
    xi_0 <- matrix(0, nrow = K, ncol = J_g)
    xi_1 <- matrix(0, nrow = K, ncol = J_g)
    for (k in 1:K) {
      xi_0[k,] <- colSums(pg_mat_g[class_assignments == k, , drop = FALSE])
      xi_1[k,] <- colSums(xi_1_raw[class_assignments == k, , drop = FALSE])
    }
    # precision_mat <- Matrix::kronecker(precision_mat, Diagonal(item_membership_list[g]))
    precision_mat <- Matrix::kronecker(Matrix::Diagonal(J_g), precision_mat)
    if (any(abs(xi_0) < 1e-5)){
      # cat("xi =", xi_0, "\n")
      xi_0[xi_0 < 1e-5] <- 0.01
    }
    # cat("dim(precision_mat) =", dim(precision_mat), "\n")
    # cat("diag(precision_mat) =", diag(precision_mat))
    Matrix::diag(precision_mat) <- Matrix::diag(precision_mat) + xi_0
    new_leaf_data[,indices] <- draw_mnorm(precision_mat, c(xi_1))
  }
  
  ## sample Polya-Gamma
  pg_mat <- matrix(rpg(N*J, 2.0*a_pg, auxiliary_mat - new_leaf_data[class_assignments,]), ncol = J)
  
  ## sample auxiliary variables from truncated normals
  auxiliary_mat[,] <- rtruncnorm(N*J, a = auxiliary_mat_range[['lb']], b = auxiliary_mat_range[['ub']],
                                 mean = new_leaf_data[class_assignments,], sd = 1/sqrt(pg_mat))
  
  return(list(leaf_data = new_leaf_data, pg_data = pg_mat, auxiliary_mat = auxiliary_mat))
}



#' Sample individual class assignments Z_i, i = 1, ..., N
#'@param data a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param leaf_data a K by J matrix of \eqn{logit(theta_{kj})}
#'@param a_pg a N by J matrix of hyperparameters of the generalized logistic distribution
#'@param auxiliary_mat a N by J matrix of truncated normal variables from previous iteration
#'@param class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@return a vector of length N, where the i-th element is the class assignment of
#'    individual i
sample_class_assignment <- function(data, leaf_data, a_pg, auxiliary_mat, class_probability){
  # number of classes
  K <- length(class_probability)
  # number of individuals
  N <- nrow(data)
  J <- ncol(leaf_data)
  # KxN log probabilities: the i-th column is the unnormalized log probabilities of K classes
  # for indiviudal i
  part1 <- log(class_probability) - a_pg * rowSums(leaf_data)
  leaf_data_t <- t(leaf_data)
  class_assignments <- rep(0, N)
  # for each individual
  for (i in 1:N) {
    logexp_sum <- colSums(log1p(exp((auxiliary_mat[i,] - leaf_data_t))))
    log_probs <- part1 - 2*a_pg * logexp_sum
    log_probs <- log_probs - max(log_probs)
    class_assignments[i] <- sample(1:K, 1, prob = exp_normalize(log_probs))
  }
  
  return(class_assignments)
}

