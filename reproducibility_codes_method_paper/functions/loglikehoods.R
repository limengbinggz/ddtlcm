
library(phylobase)
library(ape)
# library(data.table)
# library(ggplot2)
# library(parallel)

# Compute the loglikelihood of a given DDT --------------------------------

#' compute the marginal log likelihood of the locations
#' tree_phylo4d: tree in the phylo4d class
#' Sigma_by_group: variance of the brownian motion
#' dist_mat_old: a known pairwise MRCA matrix between leaf nodes, to save computation time
logllk_location <- function(tree_phylo4d, Sigma_by_group, num_items_per_group, dist_mat = NULL, tol = 1e-7){
  # get the number of items
  J <- sum(num_items_per_group)
  # get the number of groups of items
  G <- length(num_items_per_group)
  # cumulative index of item groups
  cum_J_index <- c(0, cumsum(num_items_per_group))
  # number of samples
  K <- nTips(tree_phylo4d)
  # extract locations
  locations <- tree_phylo4d@data[as.character(1:K), paste0("x", 1:J)]
  rownames(locations) <- tipLabels(tree_phylo4d)
  locations <- locations[order(as.numeric(substring(rownames(locations),2))),]
  root_node <- rootNode(tree_phylo4d)
  root_child <- names(phylobase::descendants(tree_phylo4d, root_node, "child"))
  # get the location of the root node
  root_location <- unlist(tree_phylo4d@data[as.character(rootNode(tree_phylo4d)), paste0("x", 1:J)])

  ### start construction of dist_mat_list
  if (is.null(dist_mat)) {
    ## now construct a pairwise distance matrix between leaf nodes
    leaf_nodes <- as.numeric(gsub("v", "", tipLabels(tree_phylo4d)))
    # create combinations of two leaf nodes including identical combinations
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
      branch_lengths <- diff(c(0, branch_lengths))
      names(branch_lengths) <- ancestors_internal

      # compute row variance of the matrix normal distribution of leaf nodes
      row_var_by_group <- unlist(t(sapply(1:length(pa_ancestors_internal),
                                          function(x) branch_lengths[x])))
      # row_var_by_group <- colSums(matrix(row_var_by_group, ncol = G))
      row_var_by_group <- sum(row_var_by_group)
      return(row_var_by_group)
    }

    unique_mrca_nodes <- unique(mrca_nodes)
    mrca_var_by_group <- lapply(setdiff(tree_phylo4d@label, "u1"), calculate_row_var_by_group)
    mrca_var_by_group <- matrix(unlist(mrca_var_by_group), ncol=1)
    # mrca_var_by_group <- matrix(unlist(mrca_var_by_group), ncol = G, nrow = 2*K-1, byrow = TRUE)
    rownames(mrca_var_by_group) <- setdiff(tree_phylo4d@label, "u1")
    row_var_by_group <- mrca_var_by_group[mrca_nodes,,drop=F]

    # note that given all other parameters, the eta at leaf nodes are independent between groups
    # so we calculate the log likelihood of leaf nodes by group and then add them up
    create_dist_mat_by_group <- function(g){
      ## now construct a pairwise distance matrix between leaf nodes
      dist_mat_by_group <- matrix(NA, nrow = K, ncol = K)
      # create distance matrix, equivalently, row covariance matrix
      dist_mat_by_group[lower.tri(dist_mat_by_group,diag = T)] <- row_var_by_group[,g]
      dist_mat_by_group[upper.tri(dist_mat_by_group)] <- t(dist_mat_by_group)[upper.tri(dist_mat_by_group)]

      # perform SVD on the dist_mat
      s <- svd(dist_mat_by_group)
      d_inv <- s$d
      # check whether dist_mat is invertible: all eigenvalues are positive, or full rank
      mat_rank <- sum(abs(s$d) > tol)
      is_invertible <- mat_rank == K

      # if dist_mat is not invertible, use density formula for a degenerate gaussian
      if (!is_invertible){
        # we find out the zero eigenvalues
        d_inv[which(abs(s$d) > tol)] <- 1 / s$d[which(abs(s$d) > tol)]
        d_inv[which(abs(s$d) <= tol)] <- 0
        # compute D^{-1/2} U^t
        Dhalf_Ut <- diag(sqrt(d_inv)) %*% t(s$u)
        # compute log determinant
        log_det <- sum(log(s$d[abs(s$d) > tol]))
      } else if (is_invertible){ # if dist_mat is invertible, just calculate its inverse
        # compute D^{-1/2} U^t
        Dhalf_Ut <- diag(sqrt(1/d_inv)) %*% t(s$u)
        # compute log determinant of Sigma
        log_det <- sum(log(s$d))
      }
      return(list(dist_mat_g = dist_mat_by_group, Dhalf_Ut = Dhalf_Ut, log_det = log_det, mat_rank = mat_rank))
    }
    dist_mat <- create_dist_mat_by_group(1)
  }
  ### end construction of dist_mat

  # compute loglikelihood
  logllk <- -sum(num_items_per_group * 0.5 * (K * log(Sigma_by_group * 2 * pi) + dist_mat$log_det))
  D_half_X <- dist_mat$Dhalf_Ut %*% as.matrix(locations)
  for (g in 1:G){
    D_half_X_g <- D_half_X[,(cum_J_index[g]+1):(cum_J_index[g+1]), drop=F]
    logllk <- logllk - 0.5 / Sigma_by_group[g] * sum(D_half_X_g * D_half_X_g)
  }

  if (is.na(logllk[1])){
    return(list(logllk = unname(logllk), dist_mat = dist_mat))
  } else {
    return(list(logllk = unname(logllk[1]), dist_mat = dist_mat))
  }
}



#' harmonic series
H_n <- function(k){
  if(k == 0) {
    return(0)
  }else{
    return(sum(1/1:k))
  }
}

#' compute factor in the exponent of the divergence time distribution
#' l = number of data points to the left
#' r = number of data points to the right
J_n <- function(l, r){
  sum(1/l:(r+l-1)) - H_n(r - 1)
  # H_n(r + l - 1) - H_n(l - 1) - H_n(r - 1)
}


#' compute loglikelihood of divergence times for a(t) = c/(1-t)
logllk_div_time_linear <- function(c, l, r, t){
  return( log(c) + (c*J_n(l, r)-1) * log(1-t) )
}

#' compute loglikelihood of divergence times for a(t) = c/(1-t)^2
logllk_div_time_quadratic <- function(c, l, r, t){
  return( log(c) - 2*log(1-t) + c*t/(1-t)*J_n(l, r) )
}


#' compute loglikelihood of the tree topology
logllk_tree_topology <- function(l, r){
  return(lfactorial(l-1) + lfactorial(r-1) - lfactorial(r+l-1))
}

#' compute the loglihood of the entire DDT
logllk_ddt <- function(c, is_linear, Sigma_by_group, tree_phylo4d, num_items_per_group,
                       tree_structure_old = NULL, dist_mat_old = NULL){
  if(is.null(tree_structure_old)){
    # obtain the tree topology
    tree <- extractTree(tree_phylo4d)
    # number of samples
    K <- nTips(tree_phylo4d)
    # get divergence times
    div_times <- c(nodeHeight(tree, nodeLabels(tree), from="root"), rep(1, K))
    names(div_times)[(K+1):(2*K)] <- 1:K
    # get the number of descends of each node
    num_descends <- c(sapply(descendants(tree, nodeLabels(tree), "tip"), length), rep(1, K))
    names(num_descends)[(K+1):(2*K)] <- 1:K
    # construct a dataframe containing the node information
    tree_structure_df <- data.frame(node = names(div_times), div_times = div_times, m_v = num_descends)

    # construct a data frame containing the edge information, divergence times (of descendends), and number of data
    #  points passing to the left and right
    tree_structure <- data.frame(tree@edge)
    # keep internal nodes only
    tree_structure <- tree_structure[tree_structure$ancestor != 0,]
    tree_structure <- tree_structure[tree_structure$descendant > K,]
    # add divergence time info
    tree_structure <- merge(x = tree_structure, y = tree_structure_df[, c("node", "div_times", "m_v")],
                            by.x = "descendant", by.y = "node", all.x = T)
    # get the number of data points passing through each internal node from the left side
    tree_structure2 <- merge(x = tree_structure, y = tree_structure[,c("ancestor", "m_v")],
                             by.x = "descendant", by.y = "ancestor", all.x=T, all.y=F)
    # remove the duplicated
    tree_structure2 <- tree_structure2[!duplicated(tree_structure2[,1:3]),]
    colnames(tree_structure2)[4:5] <- c("m_v","l")
    # replace NA's with 1. These nodes are the parents of leaves
    tree_structure2$l[tree_structure2$m_v==2] = 1
    # get the number of data points passing to the right
    tree_structure2$r <- tree_structure2$m_v - tree_structure2$l

    # compute the log likelihood of tree topology
    tree_structure2$logllk_tree_topology <- apply(tree_structure2, 1, function(x) logllk_tree_topology(l = x[5], r=x[6]))
    # compute the log likelihood of divergence times
    if(is_linear){
      tree_structure2$logllk_div_time <- apply(tree_structure2, 1,
                                               function(x) logllk_div_time_linear(c = c, l = x[5], r = x[6], t = x[3]))
    } else{
      tree_structure2$logllk_div_time <- apply(tree_structure2, 1,
                                               function(x) logllk_div_time_quadratic(c = c, l = x[5], r = x[6], t = x[3]))
    }
  }else{
    tree_structure2 = tree_structure_old
  }
  # compute the log likelihood of locations
  logllk_location_value <- logllk_location(tree_phylo4d =  tree_phylo4d, Sigma_by_group = Sigma_by_group,
                                           num_items_per_group, dist_mat_old)
  # cat("topology: ", sum(tree_structure2$logllk_tree_topology, tree_structure2$logllk_div_time))
  # cat("locations: ", logllk_location_value$logllk)
  logllk <- sum(tree_structure2$logllk_tree_topology, tree_structure2$logllk_div_time) + logllk_location_value$logllk
  return(list(logllk = logllk, tree_structure = tree_structure2, dist_mat = logllk_location_value$dist_mat))
}



#'@description calculate loglikelihood of the latent class model, conditional on tree structure
#'@param response_matrix a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param leaf_data a K by J matrix of logit(theta_{kj})
#'@param prior_class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@param prior_dirichlet a vector of length K. The Dirichlet prior of class probabilities
#'@param ClassItem a K by J matrix, where the k,j-th element counts the number of individuals
#'    that belong to class k have a positive response to item j
#'@param Class_count a length K vector, where the k-th element counts the number of individuals
#'    belonging to class k
#'@return a numeric of loglikelihood
logllk_lcm <- function(response_matrix, leaf_data, prior_class_probability,
                       prior_dirichlet, ClassItem, Class_count){
  # loglikelihood of the responses
  logllk <- sum(ClassItem * leaf_data) - colSums(t(log1p(exp(leaf_data)))) %*% Class_count
  # loglikelihood of the responses
  logllk <- logllk + sum((Class_count + prior_dirichlet - 1) * log(prior_class_probability))
  return(logllk[1])
}



#'@description calculate loglikelihood of the DDT-LCM
#'@param response_matrix a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param leaf_data a K by J matrix of logit(theta_{kj})
#'@param prior_class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@param prior_dirichlet a vector of length K. The Dirichlet prior of class probabilities
#'@param ClassItem a K by J matrix, where the k,j-th element counts the number of individuals
#'    that belong to class k have a positive response to item j
#'@param Class_count a length K vector, where the k-th element counts the number of individuals
#'    belonging to class k
#'@return a numeric of loglikelihood
logllk_ddt_lcm <- function(c, Sigma_by_group, tree_phylo4d, num_items_per_group,
                           tree_structure_old = NULL, dist_mat_old = NULL,
                           response_matrix, leaf_data, prior_class_probability,
                           prior_dirichlet, ClassItem, Class_count){
  logllk_ddt(c, Sigma_by_group, tree_phylo4d, num_items_per_group,
             tree_structure_old, dist_mat_old)[1] +
    logllk_lcm(response_matrix, leaf_data, prior_class_probability,
               prior_dirichlet, ClassItem, Class_count)$logllk
}




