
library(phylobase)
library(ape)
library(data.table)
library(ggplot2)
library(label.switching)

# Posterior summary functions for DDT --------------------------------

#' compute Posterior co-clustering probability
#' @param result a list of results from MCMC chain
#' @param t_vec a vector of times between 0 and 1, at which we compute PCP(t)
#' @param n_way an integer indicating how many leaves to co-cluster
#' @param mc.ncores an integer specifying the number of cores in mclapply to compute PCP over posterior samples in parallel
#' @export
#' @example
#' t_vec <- seq(0, 1, length.out = 11)
#' n_way = 3
#' mc.cores = 3L
#' pcp(result, t_vec, n_way, mc.cores)
require(parallel)
pcp <- function(result, t_vec, n_way, mc.cores = 3L){
  if (class(result) != "ddt_lcm"){
    stop("Result should be a ddt_lcm class object.")
  }
  # if the input is a phylo4d tree
  if (class(result) == "phylo4d"){
    tree_data_list = list()
    tree_data_list[[1]] <- result
    result <- list(tree_data_list = tree_data_list)
    mc.cores = 1L
  }
  # create a data frame where each row contains one of the n_way combinations of leaf nodes
  post_tree <- result$tree_samples$tree_list[[1]]
  # get divergence times of all nodes on the tree
  div_times <- nodeHeight(post_tree, post_tree@label, "root")
  # unique_div_times <- sort(unique(div_times))
  # get leaf node labels
  leaf_nodes <- grep("v", post_tree@label, value = TRUE)
  leaf_combs <- data.frame(t(combn(leaf_nodes, n_way)))
  colnames(leaf_combs) <- paste0("leaf", 1:n_way)
  # # initialize a pcp matrix, where each row represents a combination, each column counts the number of
  # # times the combination is co-clustered before time t
  # pcp_mat <- matrix(0, nrow = nrow(leaf_combs), ncol = length(t_vec))
  # colnames(pcp_mat) <- as.character(t_vec)

  calculate_pcp <- function(one_tree){
    # initialize a pcp matrix, where each row represents a combination, each column counts the number of
    # times the combination is co-clustered before time t
    pcp_mat <- matrix(0, nrow = nrow(leaf_combs), ncol = length(t_vec))
    colnames(pcp_mat) <- as.character(t_vec)
    # get the MRCA of each leaf node combination
    mrcas <- one_tree@label[apply(leaf_combs, 1, function(x) MRCA(one_tree, unlist(x)))]
    # distances to root node
    dist <- nodeHeight(one_tree, mrcas, "root")
    pcp_mat <- pcp_mat + t(sapply(dist, function(x) x > t_vec))
    return(pcp_mat)
  }
  pcp_mat <- mclapply(result$tree_samples$tree_list, calculate_pcp, mc.cores = mc.cores)
  pcp_mat <- Reduce("+", pcp_mat)


  # for (i in 1:length(result$tree_data_list)) {
  #   post_tree <- result$tree_data_list[[i]]
  #   # get the MRCA of each leaf node combination
  #   mrcas <- post_tree@label[apply(leaf_combs, 1, function(x) MRCA(post_tree, unlist(x)))]
  #   # distances to root node
  #   dist <- nodeHeight(post_tree, mrcas, "root")
  #   pcp_mat <- pcp_mat + t(sapply(dist, function(x) x < t_vec))
  # }
  # divide by the number of samples to calculate co-clustering probabilities
  pcp_mat = pcp_mat / length(result$tree_samples$tree_list)
  pcp_mat <- cbind(leaf_combs, pcp_mat)
  return(pcp_mat)
}


#' compute integrated Posterior co-clustering probability
#' @param pcp_mat a PCP matrix from the pcp() function. The first n columns indicate the n-way leaf node
#'    combinations, and the remaining columns indicate the PCP at times of the column names
#' @param mc.ncores an integer specifying the number of cores in mclapply to compute PCP over
#'    posterior samples in parallel
#' @export
#' @example
#' mc.cores = 3L
#' ipcp(pcp_mat, mc.cores)
ipcp <- function(pcp_mat, mc.cores = 3L){
  # get n-way combinations
  n_way <- length(grep("leaf", colnames(pcp_mat)))
  # get the pcp value matrix
  pcp_mat_values <- pcp_mat[-(1:n_way)]
  # get the times at which PCP were evaluated at
  t_vec <- sort(colnames(pcp_mat)[-(1:n_way)])
  # sort columns by times from small to large
  pcp_mat_values <- pcp_mat_values[, t_vec]
  # convert to numeric times
  t_vec <- as.numeric(t_vec)
  # if 1 is not in the times, then add it
  if (abs(t_vec[length(t_vec)] - 1) > 1e-6){
    t_vec <- c(t_vec, 1)
  }
  # if 0 is not in the times, then add it
  if (abs(t_vec[1]) > 1e-6){
    t_vec <- c(0, t_vec)
    pcp_mat_values <- cbind("0" = 1.0, pcp_mat_values)
  }
  # interval lengths between two times
  interval_length <- diff(t_vec)
  # if the last evaluated time is 1, then we remove it
  if (abs(t_vec[length(t_vec)] - 1) < 1e-6){
    pcp_mat_values <- pcp_mat_values[, -length(t_vec)]
  }
  ipcp_values <- apply(pcp_mat_values, 1, function(x) sum(x * interval_length))
  pcp_mat <- cbind(pcp_mat, ipcp = ipcp_values)
  return(pcp_mat)
}




#' compute integrated Posterior co-clustering probability
#' @param tree_phylo4d a phylo4d object
#' @param n_way an integer indicating how many leaves to co-cluster
#' @export
#' @example
#' mrca_distances(tree_phylo4d, n_way = 2L)
mrca_distances <- function(tree_phylo4d, n_way = 2L){
  # create a data frame where each row contains one of the n_way combinations of leaf nodes
  # get leaf node labels
  leaf_nodes <- grep("v", tree_phylo4d@label, value = TRUE)
  leaf_combs <- data.frame(t(combn(leaf_nodes, n_way)))
  colnames(leaf_combs) <- paste0("leaf", 1:n_way)

  # get the MRCA of each leaf node combination
  mrcas <- tree_phylo4d@label[apply(leaf_combs, 1, function(x) MRCA(tree_phylo4d, unlist(x)))]
  # distances to root node
  dist <- nodeHeight(tree_phylo4d, mrcas, "root")
  leaf_combs <- cbind(leaf_combs, mrca_distance = dist)

  return(leaf_combs)
}



#' compute integrated Posterior co-clustering probability
#' @param tree1 the first phylo4d object
#' @param tree2 the second phylo4d object. Note that tree1 and tree 2 should have the same number of nodes
#'     and same node labels
#' @param n_way an integer indicating how many leaves to assess MRCA branching times
#' @param difference_norm a string specifying if the difference between two upper triangular matrix should
#'     be calculated by the Frobenious (L2) norm or L0 (max) norm
#' @export
#' @example
#' mrca_distances(tree_phylo4d, n_way = 2L)
compare_mrca_distances <- function(tree1, tree2, n_way = 2L, difference_norm = c("Frobenious", "L0")){

  if ( !all(sort(tree1@label) == sort(tree1@label)) ) {
    stop("tree1 and tree2 should have equal number of nodes, and node labels should also be the same.")
  }
  difference_norm = match.arg(difference_norm)
  # compute the MRCA divergence times
  dist1 <- mrca_distances(tree1, n_way)
  dist2 <- mrca_distances(tree2, n_way)
  if (difference_norm == "Frobenious"){
    value <- sqrt(sum(dist1$mrca_distance - dist2$mrca_distance)^2)
  }else if (difference_norm == "L0"){
    value <- max(abs(dist1$mrca_distance - dist2$mrca_distance))
  }

  return(value)
}



#' plot the KxJ response probability matrix
#' @param response_probs_data a data.frame containing 3 variables:
#'     variable - categorical item index of J items
#'     value - numeric value response probability
#'     class - integer index of latent classes
#' @param num_items_per_group a vector of length G, each denoting the number of items in group g
#' @export
plot_response_probs <- function(response_probs_data,  num_items_per_group){
  p <- ggplot(data = response_probs_diff, aes(x = variable, y = class, fill = value)) +
    geom_tile(hjust = 0.5, vjust = 0.5, interpolate = FALSE, na.rm = TRUE) +
    labs(x = "Item", y = "Class", fill = "estimated - true") +
    scale_fill_gradientn(limits = c(-1,1), colours=c("gold", "white", "blue2")) +
    theme(#legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line=element_blank(),
      axis.text.y=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(size = 20),
      legend.text = element_text(size = 20),
      legend.title = element_text(size = 20)) +
    coord_cartesian(clip = "off") #+
    # theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
  cum_J_index <- c(0, cumsum(num_items_per_group))
  for (g in 1:G) {
    # print(c(cum_J_index[g]+1, cum_J_index[g+1]))
    p <- p + geom_brace(aes_(x=c(cum_J_index[g]+1, cum_J_index[g+1]), y=c(0.4, 0), label = paste("Group", g)),
                        inherit.data=F, rotate=180, labelsize=5)
  }
  return(p)
}




#' Summarize the output of a ddt_lcm model
#' @param model a ddt_lcm object
#' @method summary ddt_lcm
#' @export
summary.ddt_lcm <- function(model, burnin = 5000, relabel = T, is_linear = T, tree_map=NULL){
  if (class(model) != "ddt_lcm"){
    stop("model should be a class ddt_lcm object.")
  }
  total_iters <- length(model$loglikelihood) - 1
  G <- dim(model$Sigma_by_group_samples)[1]
  Sigma_by_group <- rep(1, G)
  num_items_per_group <- model$setting$num_items_per_group
  J <- sum(num_items_per_group)
  K <- nrow(model$class_probs_samples)

  ## get the MAP tree estimate
  map_index <- which.max(model$loglikelihood[-(1:(burnin))])
  if (length(model$tree_samples$tree_list) > 1){
    tree_map <- model$tree_samples$tree_list[[map_index+burnin]]
    # get the information of the initial tree
    tree_info <- logllk_ddt(1, is_linear, Sigma_by_group, tree_map, num_items_per_group,
                            tree_structure_old = NULL, dist_mat_old = NULL)
    tree_Sigma <- tree_info$dist_mat$dist_mat_g
  } else{
    tree_map <- NULL
    tree_Sigma <- NULL
  }
  # else{
  #   tree_map <- model$tree_samples$tree_list[[1]]
  # }

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
  Sigma_summary <- posterior_summary(x=t(model$Sigma_by_group_samples[,(burnin+1):total_iters,drop=F]), var_names)
  # divergence parameter
  c_summary <- posterior_summary(model$c_samples[(burnin+1):total_iters,drop=F], "c")

  response_probs_samples <- array(t(model$response_probs_samples[,(burnin+1):total_iters]),
                                  dim = c(total_iters-burnin, K, J))
  class_probs_samples <- t(model$class_probs_samples)[(burnin+1):total_iters,,drop=F]

  ### if we need to relabel
  if (relabel) {
    map_index <- which.max(model$loglikelihood_lcm[(burnin+1):total_iters:total_iters]) + burnin

    require(label.switching)
    # switch labels. m = # of iterations
    ls_lcm <- label.switching(
      method = c("ECR-ITERATIVE-1"),
      zpivot = model$Z_samples[,map_index],
      # mxN integer array of the latent allocation vectors generated from an MCMC algorithm
      z = t(model$Z_samples[,(burnin+1):total_iters]),
      K = K,
      # KxJ array containing the parameter that will be used as a pivot
      prapivot = matrix(model$response_probs_samples[,map_index], nrow=K),
      constraint = 1,
      mcmc = response_probs_samples
    )
    for (iter in 1:(total_iters-burnin)){
      response_probs_samples[iter,,] <- response_probs_samples[iter, ls_lcm$permutations$`ECR-`[iter,],]
      class_probs_samples[iter,] <- class_probs_samples[iter, ls_lcm$permutations$`ECR-ITERATIVE-1`[iter,]]
    }
  }
  # latent class proportions
  class_probs_summary <- posterior_summary(class_probs_samples, paste0("pi_", 1:K))
  # class response profiles
  ### need to update the extraction of J_g from model ###
  # var_names <- data.frame(k=rep(1:K, each = J), g=rep(1:G, J_g), j=unlist(sapply(J_g, seq)))
  var_names <- data.frame(k=rep(1:K, J), g=rep(1:G, num_items_per_group), j=c(unlist(sapply(num_items_per_group, seq))))
  var_names <- paste0("theta_", apply(var_names, 1, paste0, collapse = ","))
  response_probs_samples <- matrix(apply(response_probs_samples, 3, c), ncol = K*J)
  response_probs_summary <- posterior_summary(response_probs_samples, var_names)
  out <- list(tree_map = tree_map, tree_Sigma = tree_Sigma, response_probs_summary = response_probs_summary,
              class_probs_summary = class_probs_summary,
              Sigma_summary = Sigma_summary, c_summary = c_summary)
  class(out) <- "summary.ddt_lcm"
  return(out)
}




#' posterior prediction of class memberships
#' @param model a ddt_lcm object
#' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' @export
predict_class <- function(ddt_lcm_summary, data){ #, burnin = 5000
  # s <- summary.ddt_lcm(model, burnin)
  class_probability <- ddt_lcm_summary$class_probs_summary[,'Mean']
  K <- length(class_probability)
  N <- nrow(data)
  leaf_data <- logit(matrix(ddt_lcm_summary$response_probs_summary[,'Mean'], nrow = K))
  log_probs <- log(class_probability) + leaf_data %*% t(data - 1) + colSums(t(log_expit(leaf_data)))
  # log_probs <- log(class_probability) + leaf_data %*% t(data) - colSums(t(log1p(exp(leaf_data))))
  # substract the max
  log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
  predictive_probs <- apply(log_probs, 2, exp_normalize)
  class_assignments <- rep(0, N)
  # for each individual
  for (i in 1:N) {
    class_assignments[i] <- sample(1:K, 1, prob = predictive_probs[,i]) #exp_normalize(log_probs[,i])
  }
  return(list(class_assignments = class_assignments, predictive_probs = t(predictive_probs)))
}



#' @description Compute information criteria for the DDT-LCM model, including the Widely Applicable Information Criterion (WAIC),
#'
#' @param model a ddt_lcm object
#' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' @param ncores an integer specifying the number of nores to
#' @export
compute_IC <- function(model, burnin = 5000, ncores = 1L){
  require(parallel)
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


#' #' @description Compute information criteria for the DDT-LCM model, including the Widely Applicable Information Criterion (WAIC),
#' #'
#' #' @param model a ddt_lcm object
#' #' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' #' @param ncores an integer specifying the number of nores to
#' #' @export
#' compute_IC <- function(model, burnin = 5000, ncores = 1L){
#'   require(parallel)
#'   num_samples <- length(result$loglikelihood[-(1:burnin)]) - 1
#'   dat <- model$data
#'   N <- nrow(model$data)
#'   K <- model$setting$K
#'
#'   ### calculate WAIC
#'   compute_posteriorllk_matrix <- function(iter){
#'     # print(iter)
#'     c <- result$c_samples[iter]
#'     Sigma_by_group <- result$Sigma_by_group_samples[iter,]
#'     tree_phylo4d <- result$tree_samples$tree_list[iter][[1]]
#'     leaf_data <- t(matrix(result$response_probs_samples[, iter], nrow = K))
#'     # loglikelihood of tree
#'     logllk_ddt_sample <- logllk_ddt(c, Sigma_by_group, tree_phylo4d, J_g,
#'                                     tree_structure_old = NULL, dist_mat_old = NULL)$logllk
#'     # marginal loglikelihood of LCM, integrated over Z
#'     class_probability <- log(result$class_probs_samples[, iter])
#'     # logllk_lcm_sample = unlist(mclapply(1:N, function(x)
#'     #   logSumExp(class_probability + colSums(dat[x,] * log(leaf_data) + (1 - dat[x,]) * log(1 - leaf_data))),
#'     #   mc.cores = 3))
#'     logllk_lcm_sample <- rep(0, N)
#'     for (i in 1:N) {
#'       logllk_lcm_sample[i] <- logSumExp(class_probability + colSums(dat[i,] * log(leaf_data) + (1 - dat[i,]) * log(1 - leaf_data)))
#'     }
#'     return(logllk_lcm_sample)
#'   }
#'   llk_matrix <- mclapply(1:num_samples+burnin, function(x) compute_posteriorllk_matrix(x), mc.cores = ncores) #num_samples
#'   llk_matrix <- matrix(unlist(llk_matrix), nrow = N, ncol = num_samples, byrow = TRUE)
#'   WAIC_result <- WAIC(llk_matrix)
#'
#'
#'   ### calculate DIC
#'   s <- summary(model, burnin)
#'   c <- s$c_summary[,"Mean"]
#'   Sigma_by_group <- s$Sigma_summary[,"Mean"]
#'   tree_phylo4d <- s$tree_map
#'   leaf_data <- matrix(s$response_probs_summary[,'Mean'], nrow = K)
#'   logllk_ddt_sample <- logllk_ddt(c, Sigma_by_group, tree_phylo4d, J_g,
#'                                   tree_structure_old = NULL, dist_mat_old = NULL)$logllk
#'   class_probability <- s$class_probs_summary[,"Mean"]
#'   logllk_lcm_sample <- 0
#'   for (i in 1:N) {
#'     logllk_lcm_sample <- logllk_lcm_sample + logSumExp(class_probability + colSums(dat[i,] * log(leaf_data) + (1 - dat[i,]) * log(1 - leaf_data)))
#'   }
#'   DIC <- 2 * (logllk_ddt_sample + logllk_lcm_sample - mean(apply(llk_matrix, 2, logSumExp)))
#'
#'   return(list(WAIC_result = WAIC_result, DIC = DIC))
#' }



#' @description Compute the Widely Applicable Information Criterion (WAIC), also known
#'   as the Widely Available Information Criterion or the Watanable-Akaike, of Watanabe (2010).
#' @param model a N x S matrix, where N is the number of individuals and S is the number of posterior samples
#' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' @param ncores an integer specifying the number of nores to
#' @export
WAIC <- function(llk_matrix) {
  require(matrixStats)
  #  log pointwise predictive density
  S <- ncol(llk_matrix)
  part <- apply(llk_matrix, 1, logSumExp)
  lppd <- - nrow(llk_matrix) * log(S) + sum(part)
  pWAIC1 <- 2 * (lppd - sum(llk_matrix) / S)
  # pWAIC2 <- sum(.rowVars(llk_matrix))
  # lppd <- sum (log(rowMeans(exp(llk_matrix))))
  # pWAIC1 <- 2*sum( log(rowMeans(exp(llk_matrix))) - rowMeans(llk_matrix) )
  pWAIC2 <- sum(rowVars(llk_matrix))
  elppd_WAIC1 <- -2*lppd + 2*pWAIC1
  elppd_WAIC2 <- -2*lppd + 2*pWAIC2
  return(list(WAIC1=elppd_WAIC1, WAIC2=elppd_WAIC2, lppd=lppd, pWAIC2=pWAIC2, pWAIC1=pWAIC1))
}




#' #' @description Compute the DIC.
#' #' @param model a N x S matrix, where N is the number of individuals and S is the number of posterior samples
#' #' @param burnin an integer specifying the number of burn-in iterations from MCMC chain
#' #' @param ncores an integer specifying the number of nores to
#' #' @export
#' DIC <- function(llk_matrix) {
#'   require(matrixStats)
#'   #  log pointwise predictive density
#'   S <- ncol(llk_matrix)
#'   part <- apply(llk_matrix, 1, logSumExp)
#'   lppd <- - nrow(llk_matrix) * log(S) + sum(part)
#'   pWAIC1 <- 2 * (lppd - sum(llk_matrix) / S)
#'   # pWAIC2 <- sum(.rowVars(llk_matrix))
#'   # lppd <- sum (log(rowMeans(exp(llk_matrix))))
#'   # pWAIC1 <- 2*sum( log(rowMeans(exp(llk_matrix))) - rowMeans(llk_matrix) )
#'   pWAIC2 <- sum(rowVars(llk_matrix))
#'   elppd_WAIC <- -2*lppd + 2*pWAIC2
#'   return(list(WAIC=elppd_WAIC, lppd=lppd, pWAIC2=pWAIC2, pWAIC1=pWAIC1))
#' }




create_leaf_cor_matrix <- function(tree_phylo4d, mc.cores = 1L){
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
  mrca_var_by_group <- unlist(mclapply(tree_phylo4d@label, calculate_row_var_by_group, mc.cores = mc.cores))
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


