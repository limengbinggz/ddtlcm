#' ##########################################################
#' Implement a Gibbs sampler for K-class latent class models with
#' independent normal priors on the logistic-transformed
#' item responsibilities. Group-specific variance parameters.
#' Z_i ~ Cat(\pi)
#' P(Y_ij = 1 \mid Z_i = k) = \theta_{kj}
#' logistic(\theta_{kj}) ~ N(0, \sigma_g^2), for item j \in group g
#' ##########################################################


library(phylobase)
library(ape)
library(data.table)
library(BayesLogit)
library(MASS)
library(extraDistr)
library(truncnorm)
library(Matrix)

# Metropolis-Hastings algorithm for DDT -----------------------------------

#' randomly detach a subtree from a given tree
random_detach_subtree <- function(tree_phylo4){
  # number of samples
  K <- nTips(tree_phylo4)
  # get the root node
  root_name <- names(rootNode(tree_phylo4))
  # child node of root
  root_child <- names(phylobase::descendants(tree_phylo4, root_name, "children"))
  # randomly sample a node that is not the root or ch(root) to detach a subtree
  detach_node <- sample( setdiff(c(paste0("v", 1:K), paste0("u", 1:K)), # tree_phylo4@label,
                                 c(root_name, root_child )),
                         1, F)
  # cat("detach_node =", detach_node)
  # the root of the detached subtree is pa(detach_node)!
  # find the parent of the detached node
  pa_detach_node <- ancestor(tree_phylo4, detach_node)
  # get edge length between the detached node and its parent
  # detach_edge_len <- edgeLength(tree_phylo4)[getEdge(tree_phylo4, detach_node, "descendant")]
  
  # divergence time of the parent node
  pa_div_time <- nodeHeight(tree_phylo4, pa_detach_node, "root")
  # divergence time of the detached node
  detach_div_time <- nodeHeight(tree_phylo4, detach_node, "root")
  # get the leaf of the detached subtree
  subtree_leaf <- names(phylobase::descendants(tree_phylo4, detach_node, "tips"))
  
  # construct the subtree
  if(length(subtree_leaf) == 1){ # if a leaf node was detached, the subtree contains a single node
    subtree <- detach_node
  }else{
    subtree <- as(subset(tree_phylo4, node.subtree = detach_node), "phylo")
  }
  # remaining tree contains the singleton from the original root
  tree_kept_leaf <- setdiff(tipLabels(tree_phylo4), subtree_leaf)
  # construct the remaining tree
  if(length(tree_kept_leaf) == 1){ # if only one leaf is kept in the tree
    tree_kept <- list(edge = matrix(c(2,1),1,2),
                      node.label = paste0("u", 1),
                      tip.label = tree_kept_leaf,
                      edge.length = 1,
                      Nnode = 1)
    class(tree_kept) <-"phylo"
  }else{
    tree_kept <- add_root(tree_old = as(subset(tree_phylo4, tips.include = tree_kept_leaf), "phylo"),
                          root_edge_length = nodeHeight(tree_phylo4, names(phylobase::MRCA(tree_phylo4, tree_kept_leaf)), from="root"),
                          root_label = paste0("u", 1),
                          leaf_label = names(rootNode(tree_phylo4)) )
  }
  return(list(tree_detached = subtree, tree_kept = tree_kept,
              pa_detach_node_label = names(pa_detach_node), pa_div_time = pa_div_time,
              detach_div_time = detach_div_time, detach_node_label = detach_node))
}

#' add a leaf branch to an existing tree tree_old
#' div_t: divergence time of the new branch
#' new_leaf_label: the label of the newly added leaf
#' where: node name of to which node in the existing tree the new leaf branch should connect to
#' position: the numerical location of the left side of the added branch
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


add_root <- function(tree_old, root_edge_length, root_label, leaf_label){
  leaf_branch <- list(edge = matrix(c(2,1), nrow = 1),
                      node.label = root_label,
                      tip.label = leaf_label,
                      edge.length = root_edge_length,
                      Nnode = 1)
  class(leaf_branch) <- "phylo"
  tree_new <- bind.tree(leaf_branch, tree_old, where=1)
  return(tree_new)
}


#' choose a divergence time, location, and tree topology to reattach a subtree
#' tree_kept: the tree to be attached to
reattach_point <- function(tree_kept, c, is_linear=T, theta=0.0, alpha=0.0){
  # plot(tree_kept, show.node.label = TRUE)
  
  # number of leaves of the remaining tree
  tree_kept_ntip <- phylobase::nTips(tree_kept)
  new_tip <- list()
  if(tree_kept_ntip == 1){
    t <- div_time(t_u = 0, m_v = 1, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
    return(list(div_time = t, root_node=tree_kept$node.label, root_child = tree_kept$tip.label, div_dist_to_root_child = 1-t))
  }else{
    tree_kept_phylo4 <- phylo4(tree_kept, check.node.labels = "keep")
    cur_n <- tree_kept_ntip
    # root node starts from time 0
    start_time <- 0
    # get root node
    root_node <- names(rootNode(tree_kept_phylo4))
    # get the index of the child of the root node
    root_child <- names(phylobase::descendants(tree_kept_phylo4, root_node, "child"))
    # get the branch length between the root and its child
    root_child_branch_len <- edgeLength(tree_kept_phylo4)[getEdge(tree_kept_phylo4, root_child)]
    dist_to_root_child <- start_time + root_child_branch_len
    while(TRUE){
      # sample a new divergence time
      div_t <- div_time(t_u = start_time, m_v = cur_n, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
      # cat("start_time = ", start_time, "\n")
      # cat("m_v = ", cur_n, "\n")
      # cat("div_t = ", div_t, "\n")
      # cat("dist_to_root_child = ", dist_to_root_child, "\n")
      # if the new divergence time happens before the root child, then we are done by diverging at div_t
      if (div_t < dist_to_root_child){
        return(list(div_time = div_t, root_node = root_node, root_child = root_child,
                    div_dist_to_root_child = dist_to_root_child - div_t))
      }
      # otherwise, we follow one of the existing paths, with probability proportional to the number of data
      # points that already traversed the path
      branch_point <- names(descendants(tree_kept_phylo4, root_child, type="child"))
      # get the number of branches from this branch point
      K <- length(branch_point)
      # get the number of samples which previously took each branch
      n_k <- unlist(lapply(phylobase::descendants(tree_kept_phylo4, branch_point, type="tip"), length))
      names(n_k) <- branch_point
      # compute the selection probabilities
      probs <- c(n_k - alpha, theta + alpha*K) / (sum(n_k) + theta)
      # select which branch to follow, or to create a new branch
      path_idx_sel <- which.max(rmultinom(1, 1, probs))
      
      if(path_idx_sel == length(probs)){
        # if we select to create a new branch, then put divergence at the root child
        return(list(div_time = dist_to_root_child, root_node = root_node, root_child = root_child,
                    div_dist_to_root_child = dist_to_root_child - div_t))
      } else{
        # if we choose one of the existing branchs, then we need to sample the divergence time again
        selected_node <- branch_point[path_idx_sel]
        # get the number of data points which previously travel the path
        m_v <- n_k[path_idx_sel]
        if(m_v == 1){
          # sample the divergence time again
          div_t <- div_time(dist_to_root_child, m_v, c, is_linear = is_linear, theta, alpha)
          return(list(div_time = div_t, root_node = root_child, root_child = selected_node,
                      div_dist_to_root_child = 1 - div_t))
        }else{
          # if more than one data point has  traversed along this path, then we need to repeat the above
          # process again: starting from the branch point, follow the sampled path, until reaching a path that
          # has been traveled by only one data point
          start_time <- dist_to_root_child
          cur_n <- m_v
          root_node <- root_child
          root_child <- selected_node
          edge_len <- edgeLength(tree_kept_phylo4)[getEdge(tree_kept_phylo4, root_child)]
          dist_to_root_child <- start_time + edge_len
        }
      }
    }
  }
}



#' randomly attach a subtree to a given DDT
#' subtree: subtree to attach to tree_kept
#' tree_kept: the tree to be attached to
#' detach_div_time: divergence time of subtree when it was extracted from the original tree
#' pa_detach_node_label: label of the parent node of the detached node
attach_subtree <- function(subtree, tree_kept, detach_div_time, pa_detach_node_label, c,
                           is_linear = T, theta=0.0, alpha=0.0){
  # plot(subtree, show.node.label = TRUE)
  # if the subtree contains a single node
  if(length(subtree) == 1){
    # print("this1")
    # select a point on the tree to attach subtree to
    new_attach_point <- reattach_point(tree_kept, c, is_linear, theta, alpha)
    tree_kept_phylo4 <- phylo4(tree_kept, check.node.labels="keep")
    # tree_kept_phylo42 <- phylo4(tree_kept, check.node.labels="drop")
    # if the new attach time is after than the detach time in the original, then resample the attach time until it is before
    while(new_attach_point$div_time >= detach_div_time){
      new_attach_point <- reattach_point(tree_kept, c=c, is_linear=is_linear, theta=theta, alpha=alpha)
    }
    # attach the subtree to the given tree as a leaf branch
    new_phylo4 <- phylo4( add_leaf_branch(tree_kept, div_t = new_attach_point$div_time, new_leaf_label = subtree,
                                          where = getNode(tree_kept_phylo4, new_attach_point$root_child),
                                          position = nodeHeight(tree_kept_phylo4, new_attach_point$root_child, "root") - new_attach_point$div_time) )
    # add label to the new attach point
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))] <- pa_detach_node_label
    
    # new_node_order <- match(new_phylo4@label, c(paste0("v", 1:K), paste0("u", 1:K)))
    
    
    # K <- nTips(new_phylo4)
    # # need to put node indices larger than the subtree one node forward
    # subtree_index <- as.numeric(gsub(".*?([0-9]+).*", "\\1", subtree))
    # which(new_phylo4@label %in% paste0("u", 2:K))#paste0("v", (subtree_index+1):K)
    # gsub("v.*?([0-9]+).*", "\\1", new_phylo4@label)
    #
    # tip_index <- which(as.numeric(names(new_phylo4@label)) <= K)
    # new_phylo4@label[tip_index] <- paste0("v", names(new_phylo4@label)[tip_index])
    # node_index <- which(as.numeric(names(new_phylo4@label)) > K)
    # new_phylo4@label[node_index] <- paste0("u", as.numeric(names(new_phylo4@label))[node_index] - K)
    # plot(new_phylo4, show.node=T)
    # new_phylo4@label[grep("v", new_phylo4@label)] <- paste0("v", 1:nTips(new_phylo4))
    return(list(new_phylo4 = new_phylo4,
                attach_root = new_attach_point$root_node,
                attach_to = new_attach_point$root_child,
                new_div_time = new_attach_point$div_time))
  }else if( phylobase::nTips(tree_kept) == 1 ){ # if the subtree has only one leaf
    # print("this2")
    # sample a new divergence time
    div_t <- div_time(t_u = 0, m_v = 1, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
    # if the new attach time is after than the detach time in the original, then resample the attach time until it is before
    while(div_t >= detach_div_time){
      div_t <- div_time(t_u = 0, m_v = 1, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
    }
    # add a root node to the subtree
    new_subtree <- add_root(tree_old = subtree,
                            root_edge_length = detach_div_time - div_t,
                            root_label = pa_detach_node_label,
                            leaf_label = names(rootNode(phylo4(subtree))))
    # add the subtree to the given tree
    new_phylo4 <- phylo4(bind.tree(tree_kept, new_subtree, where = 1, position = 1 - div_t))
    # add label to the new attach point
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))] <- pa_detach_node_label
    return(list(new_phylo4 = new_phylo4,
                attach_root = tree_kept$node.label,
                attach_to = tree_kept$tip.label,
                new_div_time = div_t))
  } else{
    # print("this3")
    # select a point on the tree to attach subtree to
    new_attach_point <- reattach_point(tree_kept, c=c, is_linear=is_linear, theta=theta, alpha=alpha)
    # print(new_attach_point)
    tree_kept_phylo4 <- phylo4(tree_kept, check.node.labels="keep")
    while (new_attach_point$div_time >= detach_div_time){
      new_attach_point <- reattach_point(tree_kept, c=c, is_linear=is_linear, theta=theta, alpha=alpha)
    }
    # add a root node to the subtree
    new_subtree <- add_root(tree_old = subtree,
                            root_edge_length = detach_div_time - new_attach_point$div_time,
                            root_label = pa_detach_node_label,
                            leaf_label = names(rootNode(phylo4(subtree))))
    # plot(new_subtree, show.node.label = TRUE)
    new_phylo4 <- phylo4(bind.tree(tree_kept, new_subtree,
                                   where = getNode(tree_kept_phylo4, new_attach_point$root_child),
                                   position = new_attach_point$div_dist_to_root_child))
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))] <- pa_detach_node_label
    # plot(new_phylo4, show.node.label = TRUE)
    return(list(new_phylo4 = new_phylo4,
                attach_root = new_attach_point$root_node,
                attach_to = new_attach_point$root_child,
                new_div_time = new_attach_point$div_time))
  }
}





proposal_log_prob <- function(old_tree_phylo4, tree_kept, old_detach_pa_div_time, old_pa_detach_node_label, old_detach_node_label,
                              new_div_time, new_attach_root, new_attach_to, c, is_linear=T){
  if (is_linear) {
    A <- A_t_linear
    a_t <- a_t_linear
  } else{
    A <- A_t_quadratic
    a_t <- a_t_quadratic
  }
  # if the given only has one leaf
  if(nTips(tree_kept)==1){
    # we only need to consider the divergence time
    q_new <- -A(c=c, new_div_time) + log(a_t(c=c, new_div_time))
    q_old <- -A(c=c, old_detach_pa_div_time) + log(a_t(c=c, old_detach_pa_div_time))
    return(c(q_new = q_new, q_old = q_old))
  }else{
    tree_kept_phylo4 <- phylo4(tree_kept)
    # get root node name
    tree_kept_root_name <- names(rootNode(tree_kept_phylo4))
    
    ## for the old tree
    # root node name of the detached tree
    old_detach_root_name <- names(ancestor(old_tree_phylo4, old_pa_detach_node_label))
    old_detach_to <- setdiff(names(phylobase::descendants(old_tree_phylo4, old_pa_detach_node_label, "children")), old_detach_node_label)
    # get the path from root to leaf on the branch that the subtree was detached
    old_path_root <- unique(c(tree_kept_root_name, names(shortestPath(tree_kept_phylo4, old_detach_root_name, tree_kept_root_name)),
                              old_detach_root_name, old_detach_to))
    if(length(old_path_root) == 2){
      m_v <- length(phylobase::descendants(tree_kept_phylo4, old_path_root[-1], "tips"))
      names(m_v) <- old_path_root[-1]
      div_t <- old_detach_pa_div_time
      names(div_t) <- "old_detach_pa_time"
    }else{
      # get the number of data points that travel through each internal and leaf node
      m_v <- sapply(phylobase::descendants(tree_kept_phylo4, old_path_root[-1], "tips"), length)
      names(m_v) <- old_path_root[-1]
      # get the divergence time of each internal node, and the detached time
      div_t <- c(nodeHeight(tree_kept_phylo4, names(m_v)[-length(m_v)], "root"), old_detach_pa_div_time)
      names(div_t) <- c(names(m_v)[-length(m_v)], "old_detach_pa_time")
    }
    # branch [uv]
    A_div_t <- A(c=c, div_t)
    A_tu <- c(0, A_div_t[-length(A_div_t)])
    A_tv <- A_div_t
    final_A_t <- sum((A_tu - A_tv) / m_v)
    frac <- sum(log(m_v[-1] / m_v[-length(m_v)]))
    final_a_t <- log( a_t(c=c, old_detach_pa_div_time) / m_v[length(m_v)] )
    q_old = sum(final_A_t, frac, final_a_t)
    
    ## for the new tree
    new_path_root <- unique(c(tree_kept_root_name, names(shortestPath(tree_kept_phylo4, new_attach_root, tree_kept_root_name)),
                              new_attach_root, new_attach_to))
    if(length(new_path_root) == 2){
      m_v <- length(phylobase::descendants(tree_kept_phylo4, new_path_root[-1], "tips"))
      names(m_v) <- new_path_root[-1]
      div_t <- new_div_time
      names(div_t) <- "new_div_time"
    }else{
      m_v <- sapply(phylobase::descendants(tree_kept_phylo4, new_path_root[-1], "tips"), length)
      names(m_v) <- new_path_root[-1]
      div_t <- c(nodeHeight(tree_kept_phylo4, names(m_v)[-length(m_v)], "root"), new_div_time)
      names(div_t) <- c(names(m_v)[-length(m_v)], "new_div_time")
    }
    A_div_t <- A(c=c, div_t)
    A_tu <- c(0, A_div_t[-length(A_div_t)])
    A_tv <- A_div_t
    final_A_t <- sum((A_tu - A_tv)/m_v)
    frac <- sum(log( m_v[-1] / m_v[-length(m_v)] ))
    final_a_t <- log( a_t(c=c, new_div_time) / m_v[length(m_v)] )
    q_new = sum(final_A_t, frac, final_a_t)
    
    return(c(q_new = q_new, q_old = q_old))
  }
}


#'@description Sample divergence function parameter c for a(t) = c / (1-t) through Gibbs sampler
#'@param shape0 shape of the Inverse-Gamma prior
#'@param rate0 rate of the Inverse-Gamma prior
#'@param tree_structure a data.frame containing the divergence times and number of
#'    data points to the left and right branches of internal nodes on the tree
#'@return a numeric value of the newly sampled c
sample_c_linear <- function(shape0, rate0, tree_structure){
  shape <- shape0 + nrow(tree_structure)
  rate <- rate0 - sum(( sapply(tree_structure$m_v-1, function(x) H_n(x)) -
                          sapply(tree_structure$l-1,function(x) H_n(x)) -
                          sapply(tree_structure$r-1,function(x) H_n(x)) ) *
                        (log(1 - tree_structure$div_times)))
  return(rgamma(1, shape = shape, rate = rate))
}

#'@description Sample divergence function parameter c for a(t) = c / (1-t)^2 through Gibbs sampler
#'@param shape0 shape of the Inverse-Gamma prior
#'@param rate0 rate of the Inverse-Gamma prior
#'@param tree_structure a data.frame containing the divergence times and number of
#'    data points to the left and right branches of internal nodes on the tree
#'@return a numeric value of the newly sampled c
sample_c_quadratic <- function(shape0, rate0, tree_structure){
  shape <- shape0 + nrow(tree_structure)
  rate <- rate0 - sum(( sapply(tree_structure$m_v-1, function(x) H_n(x)) -
                          sapply(tree_structure$l-1,function(x) H_n(x)) -
                          sapply(tree_structure$r-1,function(x) H_n(x)) ) *
                        tree_structure$div_times / (1 - tree_structure$div_times) )
  return(rgamma(1, shape = shape, rate = rate))
}



#'@description Sample item group-specific variances through Gibbs sampler
#'@param shape0 a vector of G elements, each being the shape of the
#'    Inverse-Gamma prior of group g
#'@param rate0 a vector of G elements, each being the rate of the
#'    Inverse-Gamma prior of group g
#'@param dist_mat a list, containing the KxK tree-structured matrix of leaf nodes,
#'    where K is the number of leaves / latent classes, and SVD components
#'@param num_items_per_group a vector of G elements, each indicating the number of
#'    items in this group
#'@param locations a KxJ matrix of leaf parameters
#'@return a numeric vector of G elements, each being the newly sampled variance
#'    of the latent location of this group
sample_sigmasq <- function(shape0, rate0, num_items_per_group, locations){
  # number of item groups
  G <- length(num_items_per_group)
  cum_J_index <- c(0, cumsum(num_items_per_group))
  # number of leaves
  K <- nrow(locations)
  Sigma_by_group <- rep(0, G)
  for (g in 1:G) {
    shape <- shape0[g] + K * num_items_per_group[g] * 0.5
    # get the location of the root node
    locations_g_diff <- locations[,(cum_J_index[g]+1):(cum_J_index[g+1]), drop=F]
    
    # trace of A^t B A = vec(B^t A)^t vec(A)
    # crossprod(c(t(inv_dist_mat) %*% X_minus_mu_g), c(X_minus_mu_g))[1]/2
    rate <- rate0[g] + 0.5 * sum( locations_g_diff **2)
    Sigma_by_group[g] <- 1.0 / rgamma(1, shape = shape, rate = rate)
  }
  return(Sigma_by_group)
}




#'@description Sample a new tree topology using Metropolis-Hastings through randomly
#'    detaching and re-attaching subtrees
#'@param tree_phylo4d_old a phylo4d object of tree from the previous iteration
#'@param Sigma_by_group a vector of diffusion variances of G groups from the previous iteration
#'@param num_items_per_group a vector of G elements, each indicating the number of
#'    items in this group
#'@param c a number of the divergence hyperparameter from the previous iteration
#'@param tree_structure_old a data.frame of tree structure from the previous iteration. Each row
#'    contains information of an internal node, including divergence times, number of data points
#'    traveling through the left and right branches
#'@param dist_mat_old a list of leaf covariance matrix from the previous iteration. The list has length
#'    G, the number of item groups
#'@return a numeric vector of G elements, each being the newly sampled variance
#'    of the latent location of this group
sample_tree_topology <- function(tree_phylo4d_old, Sigma_by_group, num_items_per_group, c, is_linear=T,
                                 tree_structure_old = NULL, dist_mat_old = NULL){
  # number of leaves / classes
  K <- nTips(tree_phylo4d_old)
  old_tree <- extractTree(tree_phylo4d_old)
  # randomly detach a subtree
  random_detach <- random_detach_subtree(tree_phylo4 = old_tree)
  # randomly attach the detached subtree
  proposal_newtree <- attach_subtree(subtree = random_detach$tree_detached, tree_kept = random_detach$tree_kept,
                                     detach_div_time = random_detach$detach_div_time,
                                     pa_detach_node_label = random_detach$pa_detach_node_label,
                                     is_linear = is_linear, c = c, theta=0, alpha = 0)
  # while(nTips(proposal_newtree$new_phylo4) != K){
  # need to avoid extremely short branch in the new proposal tree to avoid degeneracy
  # in the covariance matrix
  while(proposal_newtree$new_div_time > 0.95){
    proposal_newtree <- attach_subtree(subtree = random_detach$tree_detached, tree_kept = random_detach$tree_kept,
                                       detach_div_time = random_detach$detach_div_time,
                                       pa_detach_node_label = random_detach$pa_detach_node_label,
                                       is_linear = is_linear, c = c, theta=0, alpha = 0)
  }
  # compute log likelihood of the old and proposal tree divergence times
  q_prob <- proposal_log_prob(old_tree_phylo4 = old_tree, tree_kept = random_detach$tree_kept,
                              old_detach_pa_div_time = random_detach$pa_div_time,
                              old_pa_detach_node_label = random_detach$pa_detach_node_label,
                              old_detach_node_label = random_detach$detach_node_label,
                              new_div_time = proposal_newtree$new_div_time,
                              new_attach_root = proposal_newtree$attach_root,
                              new_attach_to = proposal_newtree$attach_to,
                              c = c, is_linear = is_linear)
  # convert to phylo4d object
  new_node_order <- match(proposal_newtree$new_phylo4@label, tree_phylo4d_old@label) # c(paste0("v", 1:K), paste0("u", 1:K))
  # reorder rows of data to from leaf to internal nodes
  tree_phylo4d_old@data <- tree_phylo4d_old@data[as.character(1:(2*K)),]
  data_new_order <- tree_phylo4d_old@data[as.character(new_node_order),]
  rownames(data_new_order) <- 1:(2*K)
  new_phylo4d <- phylo4d(proposal_newtree$new_phylo4, all.data=data_new_order)#tree_phylo4d_old@data
  # new_phylo4d <- phylo4d(proposal_newtree$new_phylo4, all.data=tree_phylo4d_old@data)#tree_phylo4d_old@data
  
  # compute log likelihood of the old and proposal tree locatioins
  old_tree_info <- logllk_ddt(c = c, is_linear, Sigma_by_group = Sigma_by_group, tree_phylo4d = tree_phylo4d_old,
                              num_items_per_group = num_items_per_group,
                              tree_structure_old = tree_structure_old, dist_mat_old = dist_mat_old)
  new_tree_info <- logllk_ddt(c = c, is_linear, Sigma_by_group = Sigma_by_group, tree_phylo4d = new_phylo4d,
                              num_items_per_group = num_items_per_group)
  
  # MH ratio
  r <- new_tree_info$logllk + q_prob["q_old"] - old_tree_info$logllk - q_prob["q_new"]
  if(log(runif(1)) < r){ # accept
    final <- list(
      accept = 1,
      logllk_model = new_tree_info$logllk,
      tree_phylo4d = new_phylo4d,
      dist_mat = new_tree_info$dist_mat,
      tree_structure = new_tree_info$tree_structure
    )
    
  }else{
    final <- list(
      accept = 0,
      logllk_model = old_tree_info$logllk,
      tree_phylo4d = tree_phylo4d_old,
      dist_mat = old_tree_info$dist_mat,
      tree_structure = old_tree_info$tree_structure
    )
  }
  return(final)
}



#'@description Sample the leaf locations and Polya-Gamma auxilliary variables
#'@param num_items_per_group a vector of G elements, each indicating the number of
#'    items in this group
#'@param root_location a vector of length J. Location of the root node on the tree
#'@param dist_mat_old a list of leaf covariance matrix from the previous iteration. The list
#'    has length G, the number of item groups
#'@param Sigma_by_group a vector of length G, each denoting the variance of the
#'    brownian motion
#'@param pg_mat a K by J matrix of PG variables from the previous iteration
#'@param Kappa a K by J matrix of parameters when using PG augmentation. It equals
#'    ClassItem - 0.5 * Class_count
#'@param Class_count a vector of length K, where the k-th element counts the number of
#'    individuals belonging to class k
#'@return a numeric vector of G elements, each being the newly sampled variance
#'    of the latent location of this group
sample_leaf_locations_pg <- function(K, num_items_per_group, Sigma_by_group,
                                     pg_mat, a_pg, auxiliary_mat, auxiliary_mat_range, class_assignments){
  N <- nrow(pg_mat)
  # total number of items
  J <- sum(num_items_per_group)
  # cumulative index of item groups
  cum_J_index <- c(0, cumsum(num_items_per_group))
  new_leaf_data <- matrix(NA, nrow = K, ncol = J)
  for (g in 1:length(num_items_per_group)) {
    precision_mat <- diag(1 / Sigma_by_group[g], K)
    # extract indices
    indices <- (cum_J_index[g]+1):(cum_J_index[g+1])
    pg_mat_g <- pg_mat[, indices, drop = F]
    auxiliary_mat_g <- auxiliary_mat[, indices, drop = F]
    # cat("g =",g, "\n")
    # cat("pg_mat_g:", pg_mat_g[1:5,], "\n")
    # cat("auxiliary_mat_g:", dim(auxiliary_mat_g), "\n")
    xi_1_raw <- pg_mat_g * auxiliary_mat_g
    xi_0 <- matrix(0, nrow = K, ncol = num_items_per_group[g])
    xi_1 <- matrix(0, nrow = K, ncol = num_items_per_group[g])
    for (k in 1:K) {
      xi_0[k,] <- colSums(pg_mat_g[class_assignments == k, , drop = F])
      xi_1[k,] <- colSums(xi_1_raw[class_assignments == k, , drop = F])
    }
    # precision_mat <- Matrix::kronecker(precision_mat, Diagonal(num_items_per_group[g]))
    precision_mat <- Matrix::kronecker(Diagonal(num_items_per_group[g]), precision_mat)
    if (any(abs(xi_0) < 1e-5)){
      # cat("xi =", xi_0, "\n")
      xi_0[xi_0 < 1e-5] <- 0.01
    }
    diag(precision_mat) <- diag(precision_mat) + xi_0
    new_leaf_data[,indices] <- draw_mnorm(precision_mat, c(xi_1))
  }

  ## sample Polya-Gamma
  pg_mat <- matrix(rpg(N*J, 2.0*a_pg, auxiliary_mat - new_leaf_data[class_assignments,]), ncol = J)

  ## sample auxiliary variables from truncated normals
  auxiliary_mat[,] <- rtruncnorm(N*J, a = auxiliary_mat_range[['lb']], b = auxiliary_mat_range[['ub']],
                                 mean = new_leaf_data[class_assignments,], sd = 1/sqrt(pg_mat))

  return(list(leaf_data = new_leaf_data, pg_data = pg_mat, auxiliary_mat = auxiliary_mat))
}



#'@description sample individual class assignment Z_i, i = 1, ..., N
#'@param data a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param leaf_data a K by J matrix of logit(theta_{kj})
#'@param class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@return a vector of length N, where the i-th element is the class assignment of
#'    individual i
sample_class_assignment <- function(data, leaf_data, class_probability){
  # number of classes
  K <- length(class_probability)
  # number of individuals
  N <- nrow(data)
  # KxN log probabilities: the i-th column is the unnormalized log probabilities of K classes
  # for indiviudal i
  log_probs <- log(class_probability) + leaf_data %*% t(data - 1) + colSums(t(log_expit(leaf_data)))
  # log_probs <- log(class_probability) + leaf_data %*% t(data) - colSums(t(log1p(exp(leaf_data))))
  # substract the max
  log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
  class_assignments <- rep(0, N)
  # for each individual
  for (i in 1:N) {
    class_assignments[i] <- sample(1:K, 1, prob = exp_normalize(log_probs[,i]))
  }
  return(class_assignments)
}


#'@description Main sampling function
gibbs_bayeslcm <- function(K, data, num_items_per_group,
                         initials = list(), priors = list(),
                         total_iters = 15, controls){
  
  # total number of items
  J <- sum(num_items_per_group)
  N <- nrow(data)
  if (ncol(data) != J){
    stop("Number of columns in data should be equal to total number of items.")
  }
  # cumulative index of item groups
  cum_J_index <- c(0, cumsum(num_items_per_group))
  item_labels <- paste0("x", 1:J)
  Sigma_by_group <- initials$Sigma_by_group
  class_assignments <- initials$class_assignments
  class_probability <- initials$class_probability
  leaf_data <- logit(initials$item_probability)
  locations_root_diff <- leaf_data
  
  # a K by J matrix, where the k,j-th element counts the number of individuals that belong to class k
  # have a positive response to item j
  # need to convert integer to double in order for rpg() to work
  Class_count <- c(tabulate(class_assignments, K)) * 1.0
  ClassItem <- matrix(0, nrow = K, ncol = J)
  for (k in 1:K) {
    ClassItem[k,] <- colSums(data[class_assignments == k,,drop=F])
  }
  # kappa matrix
  Kappa <- ClassItem - 0.5 * Class_count
  a_pg <- priors$a_pg# 1.0
  # initial values of polya gamma variables: NxJ
  pg_mat <- matrix(NA, nrow = N, ncol = J)
  pg_mat[,] <- rpg(N*J, 2*a_pg, 0.0)
  # initialize the auxiliary variables from truncated normals
  auxiliary_mat <- matrix(0.0, nrow = N, ncol = J)
  auxiliary_mat_range <- list()
  auxiliary_mat_range[['lb']] <- c(-Inf, 0)[data+1]
  auxiliary_mat_range[['ub']] <- c(0, Inf)[data+1]
  auxiliary_mat[,] <- rtruncnorm(N*J, a = auxiliary_mat_range[['lb']], b = auxiliary_mat_range[['ub']],
                                 mean = leaf_data[class_assignments,], sd = pg_mat)
  
  # extract prior hyperparameters
  shape_c <- priors$shape_c
  rate_c <- priors$rate_c
  shape_sigma <- priors$shape_sigma
  rate_sigma <- priors$rate_sigma
  prior_dirichlet <- priors$prior_dirichlet
  prior_selection_probability <- priors$prior_selection_probability
  
  leaf_node_labels <- paste0("v", 1:K)
  
  # initialize storage
  Sigma_by_group_samples <- matrix(NA, length(Sigma_by_group), total_iters)
  response_probs_samples <- matrix(NA, nrow = K*J, ncol = total_iters)
  class_probs_samples <- matrix(0, nrow = K, ncol = total_iters)
  Z_samples <- matrix(0, nrow = N, ncol = total_iters)
  
  # model log likelihoods
  logllk_model_lcm <- rep(0, total_iters)
  
  # set.seed(2022)
  for(iter in 1:total_iters){
    if (iter %% 100 == 1){
      cat("------ iteration ", iter-1, "------\n")
    }
    
    ### Sample group variance
    Sigma_by_group <- sample_sigmasq(shape0 = shape_sigma, rate0 = rate_sigma, 
                                     num_items_per_group, locations_root_diff)
    Sigma_by_group_samples[,iter] <- Sigma_by_group
    
    
    ### Sample leaf locations and polya-gamma auxilliary variables
    # want: leaf_data: rows are v1 to vK
    # Kappa, pg_data, Class_count: 1 to K
    # but dist_mat_old: new order
    
    sampled <- sample_leaf_locations_pg(K, num_items_per_group, Sigma_by_group,
                                        pg_mat, a_pg, auxiliary_mat, auxiliary_mat_range, class_assignments)
    # replace old leaf locations with the newly sampled ones
    leaf_data <- sampled$leaf_data
    rownames(leaf_data) <- leaf_node_labels
    colnames(leaf_data) <- item_labels
    pg_mat <- sampled$pg_data
    auxiliary_mat <- sampled$auxiliary_mat
    locations_root_diff <- leaf_data

    ### Sample individual class assignments Z_i
    # class_assignments <- sample_class_assignment(data, leaf_data, a_pg, auxiliary_mat, class_probability)
    class_assignments <- sample_class_assignment(data, leaf_data, class_probability)
    Class_count <- tabulate(class_assignments, nbins = K) #c(table(class_assignments)) * 1.0
    ClassItem <- matrix(0, nrow = K, ncol = J)
    # print(Class_count)
    for (k in 1:K) {
      class_k <- data[class_assignments == k,,drop = F]
      if (length(class_k) > 0){
        ClassItem[k,] <- colSums(data[class_assignments == k,,drop = F])
      } else{
        ClassItem[k,] <- 0
      }
    }
    
    ### sample class probability
    class_probability <- c(rdirichlet(1, prior_dirichlet + Class_count))

    response_probs_samples[,iter] <- expit(leaf_data)
    class_probs_samples[,iter] <- class_probability
    Z_samples[, iter] <- class_assignments
    
    ### compute model loglikelihood
    logllk_model_lcm[iter] <- logllk_lcm(data, leaf_data, class_probability, prior_dirichlet, ClassItem, Class_count)
  }
  
  
  # combine all results in a list
  setting <- list(K = K, num_items_per_group = num_items_per_group, G = G)
  result <- list(
    Sigma_by_group_samples = Sigma_by_group_samples,
    response_probs_samples = response_probs_samples,
    class_probs_samples = class_probs_samples,
    Z_samples = Z_samples,
    loglikelihood_lcm = logllk_model_lcm,
    setting = setting,
    controls = controls,
    data = data
  )
  class(result) <- "gibbs_bayeslcm"
  return(result)
}



#' Summarize the output of a gibbs_bayeslcm model
#' @param model a gibbs_bayeslcm object
#' @method summary gibbs_bayeslcm
#' @export
summary.gibbs_bayeslcm <- function(model, burnin = 5000, relabel = T){
  if (class(model) != "gibbs_bayeslcm"){
    stop("model should be a class gibbs_bayeslcm object.")
  }
  total_iters <- length(model$loglikelihood_lcm)
  G <- dim(model$Sigma_by_group_samples)[1]
  num_items_per_group <- model$setting$num_items_per_group
  J <- sum(num_items_per_group)
  K <- model$setting$K

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
  response_probs_samples <- array(t(model$response_probs_samples[,(burnin+1):total_iters]),
                                  dim = c(total_iters-burnin, K, J))
  class_probs_samples <- t(model$class_probs_samples)[(burnin+1):total_iters,,drop=F]
  
  ### if we need to relabel
  if (relabel) {
    # print("rel")
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
  # var_names <- data.frame(k=rep(1:K, J), g=rep(1:G, num_items_per_group), j=c(unlist(sapply(num_items_per_group, seq))))
  var_names <- data.frame(k=rep(1:K, J), g=rep(1:G, num_items_per_group), j=c(unlist(sapply(num_items_per_group, seq))))
  var_names <- paste0("theta_", apply(var_names, 1, paste0, collapse = ","))
  response_probs_samples <- matrix(apply(response_probs_samples, 3, c), ncol = K*J)
  response_probs_summary <- posterior_summary(response_probs_samples, var_names)
  out <- list(response_probs_summary = response_probs_summary,
              class_probs_summary = class_probs_summary,
              Sigma_summary = Sigma_summary)
  class(out) <- "summary.gibbs_bayeslcm"
  return(out)
}







