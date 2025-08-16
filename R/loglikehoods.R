#' Compute log likelihood of parameters
#' @description Compute the marginal log likelihood of the parameters on the leaves of a tree
#' @param tree_phylo4d a "phylo4d" object
#' @param Sigma_by_group a vector of diffusion variances of G groups 
#' @param item_membership_list a list of G elements, where the g-th element contains the column
#'  indices of `data` corresponding to items in major group g
#' @param dist_mat a tree-structured covariance matrix from a given tree. Default is NULL. If given
#'  a matrix, then computation of the covariance matrix will be skipped to save time. This is useful
#'  in the Metropolis-Hasting algorithm when the previous proposal is not accepted.
#' @param tol a small number to prevent underflow when computing eigenvalues
#' @return A list of two elements: a numeric loglikelihood, a covariance matrix of the input tree
#' @family likelihood functions
logllk_location <- function(tree_phylo4d, Sigma_by_group, item_membership_list, dist_mat = NULL, tol = 1e-7){
  # phylobase::phylobase.options(singleton="ok")
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  # # cumulative index of item groups
  # cum_J_index <- c(0, cumsum(num_items_per_group))
  # number of samples
  K <- nTips(tree_phylo4d)
  # extract locations
  locations <- tree_phylo4d@data[as.character(1:K), paste0("x", 1:J)]
  rownames(locations) <- tipLabels(tree_phylo4d)
  locations <- locations[order(as.numeric(substring(rownames(locations),2))),]
  root_node <- names(rootNode(tree_phylo4d))
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
      row_var_by_group <- unlist(t(sapply(seq_along(pa_ancestors_internal),
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
      dist_mat_by_group[lower.tri(dist_mat_by_group,diag = TRUE)] <- row_var_by_group[,g]
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
  logllk <- -sum(length(item_membership_list) * 0.5 * (K * log(Sigma_by_group * 2 * pi) + dist_mat$log_det))
  D_half_X <- dist_mat$Dhalf_Ut %*% as.matrix(locations)
  for (g in 1:G){
    D_half_X_g <- D_half_X[, item_membership_list[[g]], drop=F]
    logllk <- logllk - 0.5 / Sigma_by_group[g] * sum(D_half_X_g * D_half_X_g)
  }

  if (is.na(logllk[1])){
    return(list(logllk = unname(logllk), dist_mat = dist_mat))
  } else {
    return(list(logllk = unname(logllk[1]), dist_mat = dist_mat))
  }
}



#' Harmonic series
#' @param k a positive integer
H_n <- function(k){
  if(k == 0) {
    return(0)
  }else{
    return(sum(1/1:k))
  }
}

#' Compute factor in the exponent of the divergence time distribution
#' @param l number of data points to the left
#' @param r number of data points to the right
J_n <- function(l, r){
  sum(1/l:(r+l-1)) - H_n(r - 1)
  # H_n(r + l - 1) - H_n(l - 1) - H_n(r - 1)
}


#' Compute loglikelihood of divergence times for a(t) = c/(1-t)
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param l number of data points to the left
#' @param r number of data points to the right
#' @param t a number in the interval (0, 1) indicating the divergence time
#' @family likelihood functions
logllk_div_time_one <- function(c, l, r, t){
  return( log(c) + (c*J_n(l, r)-1) * log(1-t) )
}

#' Compute loglikelihood of divergence times for a(t) = c/(1-t)^2
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param l number of data points to the left
#' @param r number of data points to the right
#' @param t a number in the interval (0, 1) indicating the divergence time
#' @family likelihood functions
logllk_div_time_two <- function(c, l, r, t){
  return( log(c) - 2*log(1-t) + c*t/(1-t)*J_n(l, r) )
}


#' Compute loglikelihood of the tree topology
#' @param l number of data points to the left
#' @param r number of data points to the right
#' @family likelihood functions
logllk_tree_topology <- function(l, r){
  return(lfactorial(l-1) + lfactorial(r-1) - lfactorial(r+l-1))
}


#' Calculate loglikelihood of a DDT, including the tree structure and node parameters
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param c_order equals 1 if using divergence function a(t) = c / (1-t), or 2 if 
#'  \eqn{a(t) = c / (1-t)^2}. Default is 1
#' @param Sigma_by_group a vector of diffusion variances of G groups from the previous iteration
#' @param tree_phylo4d a "phylo4d" object
#' @param item_membership_list a list of G elements, where the g-th element contains the column
#'  indices of `data` corresponding to items in major group g
#' @param tree_structure_old a list of at least named elements: loglikelihoods of the input tree topology
#'  and divergence times. These can be directly obtained from the return of this function.
#'  Default is NULL. If given a list, then computation of the loglikelihoods will be skipped to save time. 
#'  This is useful in the Metropolis-Hasting algorithm when the previous proposal is not accepted.
#' @param dist_mat_old a tree-structured covariance matrix from a given tree. Default is NULL. 
#' @return a numeric of loglikelihood
#' @family likelihood functions
logllk_ddt <- function(c, c_order, Sigma_by_group, tree_phylo4d, item_membership_list,
                       tree_structure_old = NULL, dist_mat_old = NULL){
  # phylobase::phylobase.options(singleton="ok")
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
                            by.x = "descendant", by.y = "node", all.x = TRUE)
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
    if(c_order){
      tree_structure2$logllk_div_time <- apply(tree_structure2, 1,
                                               function(x) logllk_div_time_one(c = c, l = x[5], r = x[6], t = x[3]))
    } else{
      tree_structure2$logllk_div_time <- apply(tree_structure2, 1,
                                               function(x) logllk_div_time_two(c = c, l = x[5], r = x[6], t = x[3]))
    }
  }else{
    tree_structure2 = tree_structure_old
  }
  # compute the log likelihood of locations
  logllk_location_value <- logllk_location(tree_phylo4d =  tree_phylo4d, Sigma_by_group = Sigma_by_group,
                                           item_membership_list, dist_mat_old)
  # cat("topology: ", sum(tree_structure2$logllk_tree_topology, tree_structure2$logllk_div_time))
  # cat("locations: ", logllk_location_value$logllk)
  logllk <- sum(tree_structure2$logllk_tree_topology, tree_structure2$logllk_div_time) + logllk_location_value$logllk
  return(list(logllk = logllk, tree_structure = tree_structure2, dist_mat = logllk_location_value$dist_mat))
}



#' Calculate loglikelihood of the latent class model, conditional on tree structure
#'@param response_matrix a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param leaf_data a K by J matrix of \eqn{logit(theta_{kj})}
#'@param prior_class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@param prior_dirichlet a vector of length K. The Dirichlet prior of class probabilities
#'@param ClassItem a K by J matrix, where the k,j-th element counts the number of individuals
#'    that belong to class k have a positive response to item j
#'@param Class_count a length K vector, where the k-th element counts the number of individuals
#'    belonging to class k
#'@return a numeric of loglikelihood
#'@family likelihood functions
logllk_lcm <- function(response_matrix, leaf_data, prior_class_probability,
                       prior_dirichlet, ClassItem, Class_count){
  # loglikelihood of the responses
  logllk <- sum(ClassItem * leaf_data) - colSums(t(log1p(exp(leaf_data)))) %*% Class_count
  # loglikelihood of the responses
  logllk <- logllk + sum((Class_count + prior_dirichlet - 1) * log(prior_class_probability))
  return(logllk[1])
}



#' Calculate loglikelihood of the DDT-LCM
#' @param leaf_data a K by J matrix of \eqn{logit(theta_kj)}
#' @param c a positive number for the divergence hyperparameter. A larger value implies
#'  earlier divergence on the tree
#' @param Sigma_by_group a vector of diffusion variances of G groups 
#' @param tree_phylo4d a "phylo4d" object
#' @param item_membership_list a list of G elements, where the g-th element contains the column
#'  indices of `data` corresponding to items in major group g
#' @param response_matrix a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#' @param tree_structure_old a list of at least named elements: loglikelihoods of the input tree topology
#'  and divergence times. These can be directly obtained from the return of this function.
#'  Default is NULL. If given a list, then computation of the loglikelihoods will be skipped to save time. 
#'  This is useful in the Metropolis-Hasting algorithm when the previous proposal is not accepted.
#' @param dist_mat_old a tree-structured covariance matrix from a given tree. Default is NULL. 
#' @param prior_class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#' @param prior_dirichlet a vector of length K. The Dirichlet prior of class probabilities
#' @param ClassItem a K by J matrix, where the k,j-th element counts the number of individuals
#'    that belong to class k have a positive response to item j
#' @param Class_count a length K vector, where the k-th element counts the number of individuals
#'    belonging to class k
#' @return a numeric of loglikelihood
#' @family likelihood functions
logllk_ddt_lcm <- function(c, Sigma_by_group, tree_phylo4d, item_membership_list,
                           tree_structure_old = NULL, dist_mat_old = NULL,
                           response_matrix, leaf_data, prior_class_probability,
                           prior_dirichlet, ClassItem, Class_count){
  logllk_ddt(c, Sigma_by_group, tree_phylo4d, item_membership_list,
             tree_structure_old, dist_mat_old)[1] +
    logllk_lcm(response_matrix, leaf_data, prior_class_probability,
               prior_dirichlet, ClassItem, Class_count)$logllk
}




