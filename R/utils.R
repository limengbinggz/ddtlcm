#' Compute divergence function
#' @description Compute value, cumulative hazard, and inverse for divergence function a(t) = c / (1-t)
#' @describeIn a_t_one value of the divergence function
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param t a number in the interval (0, 1) indicating the divergence time
#' @family divergence functions
#' @return The value and cumulative hazard return a positive number. The inverse function returns a number in the interval (0, 1).
#' @examples
#' a_t_one(1, 0.5)
#' @export
a_t_one <- function(c, t) c/(1-t)

#' @describeIn a_t_one cumulative hazard function
#' @examples
#' a_t_one_cum(1, 0.5)
#' @export
a_t_one_cum <- function(c, t) -c * log(1-t)

#' @param y a positive number to take inverse 
#' @describeIn a_t_one inverse function
#' @examples
#' A_t_inv_one(1, 2)
#' @export
A_t_inv_one <- function(c, y) 1.0 - exp(- y/c)

#' Compute divergence function
#' @description Compute value, cumulative hazard, and inverse for divergence function a(t) = c / (1-t)^2
#' @describeIn a_t_two value of the divergence function
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param t a number in the interval (0, 1) indicating the divergence time
#' @family divergence functions
#' @return The value and cumulative hazard return a positive number. The inverse function returns a number in the interval (0, 1).
#' @examples
#' a_t_two(1, 0.5)
#' @export
a_t_two <- function(c, t) c/(1-t)^2

#' @describeIn a_t_two cumulative hazard function
#' @examples
#' a_t_two_cum(1, 0.5)
#' @export
a_t_two_cum <- function(c, t) -c + c / (1.0-t)

#' @param y a positive number to take inverse 
#' @describeIn a_t_two inverse function
#' @examples
#' A_t_inv_two(1, 2)
#' @export
A_t_inv_two <- function(c, y) y / (c+y)

#' Sample divergence time on an edge uv previously traversed by m(v) data points
#' @param t_u a number in the interval (0, 1) indicating the divergence time at node u
#' @param m_v an integer for the number of data points traversed through node v
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param c_order equals 1 if using divergence function a(t) = c / (1-t), or 2 if 
#'  a(t) = c / (1-t)^2. Default is 1
#' @param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0. For general multifurcating tree from a Pitman-Yor process,
#'    specify positive values to alpha and theta. It is, however, recommended using alpha = 
#'    theta = 0 in inference because multifurcating trees have not been tested rigorously.
#' @return a number in the interval (0, 1)
div_time <- function(t_u, m_v, c, c_order = 1, alpha = 0, theta = 0){
  u <- runif(1)
  if (c_order == 1){
    x = a_t_one_cum(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u)
    return (A_t_inv_one(c, x))
  } else if (c_order == 2){
    x = a_t_two_cum(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u);
    return (A_t_inv_two(c, x))
  } else {
    stop("c_order must take value 1 or 2.")
  }
}


#' Add a leaf branch to an existing tree tree_old
#' @param tree_old the original "phylo" tree (with K leaves) to which the leaf branch will be added
#' @param div_t divergence time of the new branch
#' @param new_leaf_label the label of the newly added leaf
#' @param where node name of to which node in the existing tree the new leaf branch should connect to
#' @param position the numerical location of the left side of the added branch
#' @return a "phylo" tree with K+1 leaves
add_leaf_branch <- function(tree_old, div_t, new_leaf_label, where, position){
  leaf_branch <- list(edge = matrix(c(2,1), nrow = 1),
                      tip.label = new_leaf_label,
                      edge.length = 1.0 - div_t, Nnode=1)
  class(leaf_branch) <- "phylo"
  tree_new <- bind.tree(tree_old, leaf_branch,
                        # where = as.numeric(substr(names(where), 2, nchar(names(where)))),
                        where = where,
                        position = position)
  return(tree_new)
}

#' Add a leaf branch to an existing tree tree_old to make a multichotomus branch
#' @param tree_old the original "phylo" tree (with K leaves) to which the leaf branch will be added
#' @param div_t divergence time of the new branch
#' @param new_leaf_label the label of the newly added leaf
#' @param where node name of to which node in the existing tree the new leaf branch should connect to
#' @return a "phylo" tree with K+1 leaves that could possibly be multichotomus
add_multichotomous_tip <- function(tree_old, div_t, new_leaf_label, where){
  branch <- list(edge = matrix(c(2,1), nrow = 1),
                 tip.label = new_leaf_label,
                 edge.length = 1 - div_t,
                 Nnode = 1)
  class(branch) <- "phylo"
  tree_new <- bind.tree(tree_old, branch, 
                        where = where
                        # where = as.numeric(substr(names(where), 2, nchar(names(where))))
                        )
  return(tree_new)
}



#' Add a singular root node to an existing nonsingular tree 
#' @param tree_old the original nonsingular "phylo" tree
#' @param root_edge_length a number in (0, 1) representing the distance 
#'  between the new and the original root nodes
#' @param root_label a character label of the new root node 
#' @param leaf_label a character label of the leaf node 
#' @return a singular "phylo" tree 
add_root <- function(tree_old, root_edge_length, root_label, leaf_label){
  leaf_branch <- list(edge = matrix(c(2,1), nrow = 1),
                      node.label = root_label,
                      tip.label = leaf_label,
                      edge.length = root_edge_length,
                      Nnode = 1)
  class(leaf_branch) <- "phylo"
  tree_new <- ape::bind.tree(leaf_branch, tree_old, where=1)
  return(tree_new)
}


#' Compute normalized probabilities: exp(x_i) / sum_j exp(x_j)
#' @param x a number or rea-valued vector
#' @return a number or rea-valued vector
exp_normalize <- function(x){
  return(exp(x - logSumExp(x)))
}


#' Efficiently sample multivariate normal using precision matrix
#'  from x ~ N(Q^{-1}a, Q^{-1}), where Q^{-1} is the precision matrix
#' @param precision_mat precision matrix Q of the multivariate normal distribution
#' @param precision_a_vec a vector a such that the mean of the multivariate normal distribution is
#'  Q^{-1}a
draw_mnorm <- function(precision_mat, precision_a_vec){
  U <- chol(precision_mat)
  b <- rnorm(nrow(precision_mat))
  backsolve(U, backsolve(U, precision_a_vec, transpose = TRUE) + b)
}


#' The expit function 
#' @description The expit function: f(x) = exp(x) / (1+exp(x)), computed
#'  in a way to avoid numerical underflow.
#' @param x a value or a numeric vector between 0 and 1 (exclusive)
#' @return a number or real-valued vector
#' @export
#' @examples
#' expit(0.2)
#' expit(c(-1, -0.3, 0.6))
expit <- function(x){
  out <- x
  positives <- which(x >= 0)
  if (length(positives) > 0){
    # print("1,")
    out[positives] <- 1 / (1 + exp(-x[positives]))
    # return(1 / (1 + exp(-x)))
  }
  if (length(positives) < length(x)) {
    # print("2,")
    negatives <- which(x < 0)
    out[negatives] <- exp(x[negatives]) / (1 + exp(x[negatives]))
    # return(exp(x) / (1 + exp(x)))
  }
  return(out)
}


#' The logistic function
#' @description The logit function: f(x) = log(x / (1/x)). Large absolute values of 
#'  x will be truncated to +/- 5 after logit transformation according to its sign.
#' @param x a value or a numeric vector between 0 and 1 (exclusive)
#' @return a number or rea-valued vector
#' @export
#' @examples
#' logit(0.2)
#' logit(c(0.2, 0.6, 0.95))
logit <- function(x){
  y <- log(x / (1-x))
  large <- which(abs(y) > 5)
  if (length(large) > 0){
    y[large] <- sign(y[large]) * 5
  }
  return(y)
}


#' Numerically accurately compute f(x) = log(x / (1/x)). 
#' @param x a value or a numeric vector between 0 and 1 (exclusive)
#' @return a number or rea-valued vector
log_expit <- function(x){
  out <- x
  idx <- which(x < -33.3)
  if (length(idx) > 0){
    out[idx] <- x[idx]
  }
  idx <- which(x >= -33 & x < -18)
  if (length(idx) > 0){
    out[idx] <- x[idx] - exp(x[idx])
  }
  idx <- which((x >= -18) & (x < 37))
  if (length(idx) > 0){
    out[idx] <- -log1p(exp(-x[idx]))
  }
  return(out)
}


#' Suppress print from cat()
#' @param x evaluation of a statement that may explicitly or implicitly involve cat()
#' @param be_quiet logical. TRUE to suppress print from cat(); FALSE to continue printing
quiet <- function(x, be_quiet=TRUE) { 
  if (be_quiet){
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }else{
    x
  }
} 



#' Create a tree-structured covariance matrix from a given tree
#' @description Retrieve the covariance matrix of leaf nodes of a DDT tree
#' @param tree_phylo4d a "phylo4d" object
#' @return a K by K covariance matrix
#' @export
#'@examples
#'# load the MAP tree structure obtained from the real HCHS/SOL data
#'data(data_synthetic)
#'# extract elements into the global environment
#'list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv()) 
#'create_leaf_cor_matrix(tree_with_parameter)
create_leaf_cor_matrix <- function(tree_phylo4d){
  K <- nTips(tree_phylo4d)
  # K <- nNodes(tree_phylo4d)
  root_node <- names(rootNode(tree_phylo4d))
  
  # create combinations of two leaf nodes including identical combinations
  # leaf_grid <- expand.grid(leaf_nodes, leaf_nodes)
  leaf_grid <- expand.grid(1:K, 1:K)
  leaf_grid <- leaf_grid[leaf_grid$Var1 >= leaf_grid$Var2, ]
  # add "v" to match leaf node labels
  leaf_grid[] <- data.frame(lapply(leaf_grid, function(x) paste0("v", x)))
  ## we will construct a row covariance matrix for the matrix normal
  # get the MRCA of each pair of leaf nodes
  mrca_nodes <- tree_phylo4d@label[
    as.character(apply(leaf_grid, 1, function(x) phylobase::MRCA(tree_phylo4d, x))) ]
  
  calculate_row_var_by_group <- function(node){
    # internal ancestors of the MRCA
    ancestors_internal <- c(names(shortestPath(tree_phylo4d, root_node, node)), node)
    # parents of the internal ancestors
    pa_ancestors_internal <- names(ancestor(tree_phylo4d, ancestors_internal))
    # get divergence times of all ancestors of the MRCA
    branch_lengths <- nodeHeight(tree_phylo4d, c(root_node, ancestors_internal), "root")
    # branch_lengths <- diff(c(0, branch_lengths))
    branch_lengths <- diff(branch_lengths)
    names(branch_lengths) <- ancestors_internal
    
    # compute row variance of the matrix normal distribution of leaf nodes
    row_var_by_group <- sum(branch_lengths)
    return(row_var_by_group)
  }
  unique_mrca_nodes <- unique(mrca_nodes)
  mrca_var_by_group <- unlist(lapply(tree_phylo4d@label, calculate_row_var_by_group))
  # mrca_var_by_group <- matrix(unlist(mrca_var_by_group), ncol = G, nrow = 2*K-1, byrow = TRUE)
  names(mrca_var_by_group) <- tree_phylo4d@label
  row_var_by_group <- mrca_var_by_group[mrca_nodes]
  
  ## now construct a pairwise distance matrix between leaf nodes
  tree_Sigma_hclust <- matrix(NA, nrow = K, ncol = K)
  # create distance matrix, equivalently, row covariance matrix
  tree_Sigma_hclust[lower.tri(tree_Sigma_hclust,diag = TRUE)] <- row_var_by_group
  tree_Sigma_hclust[upper.tri(tree_Sigma_hclust)] <- t(tree_Sigma_hclust)[upper.tri(tree_Sigma_hclust)]
  rownames(tree_Sigma_hclust) <- colnames(tree_Sigma_hclust) <- paste0("v", 1:K)
  
  return(tree_Sigma_hclust)
}




#' Compute information criteria for the DDT-LCM model
#' @description Compute information criteria for the DDT-LCM model, including the Widely Applicable 
#'  Information Criterion (WAIC), and Deviance Information Criterion (DIC). WAIC and DIC are computed 
#'  using two different methods described in Gelman, Hwang, and Vehtari (2013), one based on (1) posterior
#'  means and the other based on (2) posterior variances.
#' @param result a "ddt_lcm" object
#' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' @param ncores an integer specifying the number of cores to compute marginal posterior log-likelihood
#'  in parallel
#' @importFrom parallel mclapply
#' @return a named list of the following elements
#' \describe{
#' \item{`WAIC_result`}{a list of WAIC-related results computed using the two methods}
#' \item{`DIC1`}{DIC computed using method 1.}
#' \item{`DIC2`}{DIC computed using method 2.}
#' }
#' @export
#' @examples
#' data(result_hchs)
#' IC_result <- compute_IC(result = result_hchs, burnin = 50, ncores = 1L)
compute_IC <- function(result, burnin = 5000, ncores = 1L){
  # require(parallel)
  num_samples <- length(result$loglikelihood[-(1:burnin)]) - 1
  dat <- result$data
  N <- nrow(result$data)
  K <- result$setting$K
  
  ### calculate WAIC
  compute_posteriorllk_matrix <- function(iter){
    # print(iter)
    # leaf_data <- t(matrix(result$response_probs_samples[iter, ,], nrow = K))
    leaf_data <- t(result$response_probs_samples[iter, ,])
    # marginal loglikelihood of LCM, integrated over Z
    class_probability <- log(result$class_probs_samples[, iter])
    logllk_lcm_sample <- rep(0, N)
    for (i in 1:N) {
      logllk_lcm_sample[i] <- logSumExp(class_probability + colSums(dat[i,] * log(leaf_data) + (1 - dat[i,]) * log(1 - leaf_data)))
    }
    return(logllk_lcm_sample)
  }
  llk_matrix <- mclapply(1:num_samples+burnin, function(x) compute_posteriorllk_matrix(x), mc.cores = ncores) #num_samples
  llk_matrix <- matrix(unlist(llk_matrix), nrow = N, ncol = num_samples, byrow = TRUE)
  WAIC_result <- WAIC(llk_matrix)
  
  
  ### calculate DIC
  s <- summary(result, burnin)
  leaf_data <- t(matrix(s$response_probs_summary[,'Mean'], nrow = K))
  class_probability <- s$class_probs_summary[,"Mean"]
  logllk_lcm_sample_thetahat <- 0
  for (i in 1:N) {
    logllk_lcm_sample_thetahat <- logllk_lcm_sample_thetahat +
      logSumExp(class_probability + colSums(dat[i,] * log(leaf_data) + (1 - dat[i,]) * log(1 - leaf_data)))
  }
  logllk_allindividuals <- colSums(llk_matrix) #apply(llk_matrix, 2, logSumExp)
  p_DIC_1 <- 2 * (logllk_lcm_sample_thetahat - mean(logllk_allindividuals))
  p_DIC_2 <- 2 * var(logllk_allindividuals)
  
  DIC1 <- -2 * (logllk_lcm_sample_thetahat + p_DIC_1)
  DIC2 <- -2 * (logllk_lcm_sample_thetahat + p_DIC_2)
  
  return(list(WAIC_result = WAIC_result, DIC1 = DIC1, DIC2 = DIC2))
}


#' Compute WAIC
#' @description Compute the Widely Applicable Information Criterion (WAIC), also known
#'   as the Widely Available Information Criterion or the Watanable-Akaike, of Watanabe (2010).
#' @param llk_matrix a N x S matrix, where N is the number of individuals and S is the number of posterior samples
#' @importFrom matrixStats rowVars
#' @return a named list
WAIC <- function(llk_matrix) {
  # require(matrixStats)
  #  log pointwise predictive density
  S <- ncol(llk_matrix)
  part <- apply(llk_matrix, 1, logSumExp)
  lppd <- - nrow(llk_matrix) * log(S) + sum(part)
  pWAIC1 <- 2 * (lppd - sum(llk_matrix) / S)
  pWAIC2 <- sum(rowVars(llk_matrix))
  elppd_WAIC1 <- -2*lppd + 2*pWAIC1
  elppd_WAIC2 <- -2*lppd + 2*pWAIC2
  return(list(WAIC1=elppd_WAIC1, WAIC2=elppd_WAIC2, lppd=lppd, pWAIC2=pWAIC2, pWAIC1=pWAIC1))
}


