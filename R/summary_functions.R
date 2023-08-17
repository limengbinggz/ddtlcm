
#' Summarize the output of a ddt_lcm model
#' @param result a "ddt_lcm" object
#' @param burnin number of samples to discard from the posterior chain as burn-ins. Default is 3000.
#' @param relabel If TRUE, perform post-hoc label switching using the Equivalence Classes 
#'  Representatives (ECR) method to solve non-identifiability issue in mixture models. If FALSE,
#'  no label switching algorithm will be performed. 
#' @param be_quiet If TRUE, do not print information during summarization. If FALSE, print label 
#'  switching information and model summary.
#' @method summary ddt_lcm
#' @importFrom label.switching label.switching
#' @family ddt_lcm results
#' @export
summary.ddt_lcm <- function(result, burnin = 3000, relabel = T, be_quiet = F){
  if (!inherits(result, "ddt_lcm")){
    stop("result should be a class ddt_lcm object.")
  }
  total_iters <- length(result$loglikelihood)
  item_membership_list <- result$setting$item_membership_list
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  Sigma_by_group <- rep(1, G)
  K <- result$setting$K
  # num_items_per_group <- unlist(lapply(1:G, function(x) rep(x, length(item_membership_list[[x]]))))
  num_items_per_group <- unlist(lapply(item_membership_list, length))

  ## get the MAP tree estimate
  map_index <- which.max(result$loglikelihood[-(1:burnin)])
  if (length(result$tree_samples$tree_list) > 1){
    tree_map <- result$tree_samples$tree_list[[map_index+burnin]]
    # get the information of the initial tree
    tree_info <- logllk_ddt(1, result$controls$c_order, Sigma_by_group, tree_map, item_membership_list,
                            tree_structure_old = NULL, dist_mat_old = NULL)
    tree_Sigma <- tree_info$dist_mat$dist_mat_g
  } else{
    tree_map <- NULL
    tree_Sigma <- NULL
  }


  ## get posterior summary of Sigma and c
  posterior_summary <- function(x, var_names){
    tryCatch({
      if (is.null(dim(x))){
        out <- matrix(0, nrow = 1, ncol = 7)
        out[,1] <- mean(x, na.rm = TRUE)
        out[,2] <- sd(x, na.rm = TRUE)
        out[,3:7] <- quantile(x, probs = c(2.5, 25, 50, 75, 97.5)/100, na.rm = TRUE)
      } else{
        out <- matrix(0, nrow = ncol(x), ncol = 7)
        out[,1] <- colMeans(x, na.rm = TRUE)
        out[,2] <- colSds(x, na.rm = TRUE)
        out[,3:7] <- t(apply(x, 2, quantile, probs = c(2.5, 25, 50, 75, 97.5)/100, na.rm = TRUE))
      }
      colnames(out) <- c("Mean", "SD", "2.5%", "25%", "50%", "75%", "97.5%")
      rownames(out) <- var_names
      return(out)
    }, error = function(cond) {
      return(NULL)
    })
  }
  var_names <- paste0("Sigma_", 1:G)
  # variances
  Sigma_summary <- posterior_summary(x=t(result$Sigma_by_group_samples[,(burnin+1):total_iters,drop=F]), var_names)
  # divergence parameter
  c_summary <- posterior_summary(result$c_samples[(burnin+1):total_iters,drop=F], "c")

  response_probs_samples <- result$response_probs_samples[(burnin+1):total_iters,,]
  # response_probs_samples <- array(t(result$response_probs_samples[,(burnin+1):total_iters]),
  #                                 dim = c(total_iters-burnin, K, J))
  class_probs_samples <- t(result$class_probs_samples)[(burnin+1):total_iters,,drop=F]

  ### if we need to relabel
  if (relabel) {
    map_index <- which.max(result$loglikelihood_lcm[(burnin+1):total_iters]) + burnin

    # require(label.switching)
    # switch labels. m = # of iterations
    # sink("NUL")
    quiet(
    ls_lcm <- label.switching(
      method = c("ECR"),
      zpivot = result$Z_samples[,map_index],
      # mxN integer array of the latent allocation vectors generated from an MCMC algorithm
      z = t(result$Z_samples[,(burnin+1):total_iters]),
      K = K,
      # KxJ array containing the parameter that will be used as a pivot
      prapivot = matrix(result$response_probs_samples[map_index,,], nrow=K),
      constraint = 1,
      mcmc = response_probs_samples
    ), be_quiet
    )
    # sink()
    for (iter in 1:(total_iters-burnin)){
      response_probs_samples[iter,,] <- response_probs_samples[iter, ls_lcm$permutations$`ECR`[iter,],]
      class_probs_samples[iter,] <- class_probs_samples[iter, ls_lcm$permutations$`ECR`[iter,]]
    }
  }
  # latent class proportions
  class_probs_summary <- posterior_summary(class_probs_samples, paste0("pi_", 1:K))
  # class response profiles
  ### need to update the extraction of J_g from result ###
  # var_names <- data.frame(k=rep(1:K, each = J), g=rep(1:G, num_items_per_group), j=unlist(sapply(num_items_per_group, seq)))
  var_names <- data.frame(k=rep(1:K, J), g=rep(1:G, num_items_per_group), j = rep(c(unlist(sapply(num_items_per_group, seq))), each = K) )
  var_names <- paste0("theta_", apply(var_names, 1, paste0, collapse = ","))
  response_probs_samples <- matrix(apply(response_probs_samples, 3, c), ncol = K*J)
  response_probs_summary <- posterior_summary(response_probs_samples, var_names)
  max_llk_full <- result$loglikelihood[map_index]
  max_llk_lcm <- result$loglikelihood_lcm[map_index]
  out <- list(tree_map = tree_map, tree_Sigma = tree_Sigma, response_probs_summary = response_probs_summary,
              class_probs_summary = class_probs_summary,
              Sigma_summary = Sigma_summary, c_summary = c_summary,
              max_llk_full = max_llk_full, max_llk_lcm = max_llk_lcm, 
              setting = result$setting, data = result$data)
  class(out) <- "summary.ddt_lcm"
  
  if (!be_quiet){
    print(out)
  }
  return(out)
}



#' Print out a ddt_lcm model
#' @param model a "ddt_lcm" object
#' @method print ddt_lcm
#' @family ddt_lcm results
#' @export
print.ddt_lcm <- function(model, digits = 3L){
  # total number of items
  J <- length(unlist(model$setting$item_membership_list))
  # number of major item groups
  G <- length(model$setting$item_membership_list)
  N <- nrow(model$data)
  K <- model$setting$K
  total_iters <- length(model$loglikelihood)

  cat("\n---------------------------------------------\n")
  cat("DDT-LCM with K =", K, "latent classes run on", N, "observations and", J, "items in", G, "major groups.\n")
  cat(total_iters, "iterations of posterior samples drawn.")
  cat("\n---------------------------------------------\n")
}


#' Print out summary of a ddt_lcm model
#' @param model_summary a "summary.ddt_lcm" object
#' @method print summary.ddt_lcm
#' @family ddt_lcm results
#' @export
print.summary.ddt_lcm <- function(model_summary, digits = 3L){
  # total number of items
  J <- length(unlist(model_summary$setting$item_membership_list))
  # number of major item groups
  G <- length(model_summary$setting$item_membership_list)
  N <- nrow(model_summary$data)
  K <- model_summary$setting$K
  total_iters <- length(model_summary$loglikelihood)
  
  cat("Conditional item response (column) probabilities of positive responses, by outcome variable, for each class (row):\n\n")
  item_names <- colnames(model_summary$data)
  response_probs_mean <- matrix(model_summary$response_probs_summary[,"Mean"], nrow = K)
  rownames(response_probs_mean) <- paste0("Class", 1:K)
  colnames(response_probs_mean) <- item_names
  print(round(response_probs_mean, digits))
  cat("\n\n")
  
  cat("Estimated class probabilities: \n", round(model_summary$class_probs_summary[,"Mean"], digits), "\n\n")
  
  cat("Maximum log-likelihood of full model:", round(model_summary$max_llk_full, digits), "\n \n")
  cat("Maximum log-likelihood of LCM only:", round(model_summary$max_llk_lcm, digits), "\n \n")
}



#' Create a tree-structured covariance matrix from a given tree
#' @description Retrieve the covariance matrix of leaf nodes of a DDT tree
#' @param tree_phylo4d a "phylo4d" object
#' @return a K by K covariance matrix
#' @export
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
  tree_Sigma_hclust[lower.tri(tree_Sigma_hclust,diag = T)] <- row_var_by_group
  tree_Sigma_hclust[upper.tri(tree_Sigma_hclust)] <- t(tree_Sigma_hclust)[upper.tri(tree_Sigma_hclust)]
  rownames(tree_Sigma_hclust) <- colnames(tree_Sigma_hclust) <- paste0("v", 1:K)

  return(tree_Sigma_hclust)
}




#' Compute information criteria for the DDT-LCM model
#' @description Compute information criteria for the DDT-LCM model, including the Widely Applicable 
#'  Information Criterion (WAIC), and Deviance Information Criterion (DIC). WAIC and DIC are computed 
#'  using two different methods described in Gelman, Hwang, and Vehtari (2013), one based on (1) posterior
#'  means and the other based on (2) posterior variances.
#' @param model a "ddt_lcm" object
#' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' @param ncores an integer specifying the number of nores to
#' @importFrom parallel mclapply
#' @return a named list of the following elements
#' \describe{
#' \item{`WAIC_result`}{a list of WAIC-related results computed using the two methods}
#' \item{`DIC1`}{DIC computed using method 1.}
#' \item{`DIC2`}{DIC computed using method 2.}
#' }
#' @export
compute_IC <- function(model, burnin = 5000, ncores = 1L){
  # require(parallel)
  num_samples <- length(model$loglikelihood[-(1:burnin)]) - 1
  dat <- model$data
  N <- nrow(model$data)
  K <- model$setting$K
  
  ### calculate WAIC
  compute_posteriorllk_matrix <- function(iter){
    # print(iter)
    leaf_data <- t(matrix(model$response_probs_samples[, iter], nrow = K))
    # marginal loglikelihood of LCM, integrated over Z
    class_probability <- log(model$class_probs_samples[, iter])
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
  s <- summary(model, burnin)
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
#' @export
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

