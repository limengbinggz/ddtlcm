###############################################################
#' Metropolis-Hasting algorithm for sampling tree topology and
#' branch lengths from the DDT branching process.
###############################################################

#' @description Randomly detach a subtree from a given tree
#' @param tree_phylo4 a "phylo4" object
#' @return a list containing the following elements:
#' \describe{
#' \item{`tree_detached`}{a "phylo" tree detached from the input tree}
#' \item{`tree_kept`}{the remaining "phylo" tree after detachment}
#' \item{`pa_detach_node_label`}{a character label of the parent of the node 
#'  from which the detachment happens}
#' \item{`pa_div_time`}{a number in (0, 1) indicating the divergence time of 
#'  the parent of the detached node}
#' \item{`detach_div_time`}{a number in (0, 1) indicating the divergence time of the detached node}
#' \item{`detach_node_label`}{a character label of the parent of the detached node}
#' }
#' @family sample trees
#' @export
#' @examples
#' library(phylobase)
#' # load the MAP tree structure obtained from the real HCHS/SOL data
#' data(data_synthetic)
#' # extract elements into the global environment
#' list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv()) 
#' detachment <- random_detach_subtree(extractTree(tree_with_parameter)) 
random_detach_subtree <- function(tree_phylo4){
  # phylobase::phylobase.options(singleton="ok")
  # number of samples
  K <- nTips(tree_phylo4)
  # get the root node
  root_name <- names(phylobase::rootNode(tree_phylo4))
  # child node of root
  root_child <- names(phylobase::descendants(tree_phylo4, root_name, "children"))
  # randomly sample a node that is not the root or ch(root) to detach a subtree
  detach_node <- sample( setdiff(c(paste0("v", 1:K), paste0("u", 1:K)), 
                                 c(root_name, root_child )),
                         1, F)
  # cat("detach_node =", detach_node)
  # the root of the detached subtree is pa(detach_node)!
  # find the parent of the detached node
  pa_detach_node <- phylobase::ancestor(tree_phylo4, detach_node)

  # divergence time of the parent node
  pa_div_time <- phylobase::nodeHeight(tree_phylo4, pa_detach_node, "root")
  # divergence time of the detached node
  detach_div_time <- phylobase::nodeHeight(tree_phylo4, detach_node, "root")
  # get the leaf of the detached subtree
  subtree_leaf <- names(phylobase::descendants(tree_phylo4, detach_node, "tips"))
  
  # construct the subtree
  if(length(subtree_leaf) == 1){ # if a leaf node was detached, the subtree contains a single node
    subtree <- detach_node
  }else{
    subtree <- as(phylobase::subset(tree_phylo4, node.subtree = detach_node), "phylo")
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
    tree_kept <- add_root(tree_old = as(phylobase::subset(tree_phylo4, tips.include = tree_kept_leaf), "phylo"),
                          root_edge_length = phylobase::nodeHeight(tree_phylo4, names(phylobase::MRCA(tree_phylo4, tree_kept_leaf)), from="root"),
                          root_label = paste0("u", 1),
                          leaf_label = names(rootNode(tree_phylo4)) )
  }
  return(list(tree_detached = subtree, tree_kept = tree_kept,
              pa_detach_node_label = names(pa_detach_node), pa_div_time = pa_div_time,
              detach_div_time = detach_div_time, detach_node_label = detach_node))
}




#' Attach a subtree to a given DDT at a randomly selected location
#' @param tree_kept the tree to be attached to
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#'  \eqn{a(t) = c/(1-t) or c/(1-t)^2}.
#' @param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0. For general multifurcating tree from a Pitman-Yor process,
#'    specify positive values to alpha and theta. It is, however, recommended using alpha = 
#'    theta = 0 in inference because multifurcating trees have not been tested rigorously.
#' @return a list of the following objects:
#' \describe{
#' \item{`div_time`}{a numeric value of newly sampled divergence time. Between 0 and 1.}
#' \item{`root_node`}{a character. Label of the root node of `tree_kept`.}
#' \item{`root_child`}{a character. Label of the child node of the root of `tree_kept`.}
#' \item{`div_dist_to_root_child`}{a N-vector with integer entries from 1, ..., K. The initial values for
#'  individual class assignments.}
#' }
#' @family sample trees
reattach_point <- function(tree_kept, c, c_order=1, theta=0.0, alpha=0.0){
  # number of leaves of the remaining tree
  tree_kept_ntip <- phylobase::nTips(tree_kept)
  new_tip <- list()
  if(tree_kept_ntip == 1){
    t <- div_time(t_u = 0, m_v = 1, c = c, c_order = c_order, theta = theta, alpha = alpha)
    return(list(div_time = t, root_node=tree_kept$node.label, root_child = tree_kept$tip.label, div_dist_to_root_child = 1-t))
  }else{
    tree_kept_phylo4 <- phylobase::phylo4(tree_kept, check.node.labels = "keep")
    cur_n <- tree_kept_ntip
    # root node starts from time 0
    start_time <- 0
    # get root node
    root_node <- names(phylobase::rootNode(tree_kept_phylo4))
    # get the index of the child of the root node
    root_child <- names(phylobase::descendants(tree_kept_phylo4, root_node, "child"))
    # get the branch length between the root and its child
    root_child_branch_len <- phylobase::edgeLength(tree_kept_phylo4)[phylobase::getEdge(tree_kept_phylo4, root_child)]
    dist_to_root_child <- start_time + root_child_branch_len
    while(TRUE){
      # sample a new divergence time
      div_t <- div_time(t_u = start_time, m_v = cur_n, c = c, c_order = c_order, theta = theta, alpha = alpha)
      # if the new divergence time happens before the root child, then we are done by diverging at div_t
      if (div_t < dist_to_root_child){
        return(list(div_time = div_t, root_node = root_node, root_child = root_child,
                    div_dist_to_root_child = dist_to_root_child - div_t))
      }
      # otherwise, we follow one of the existing paths, with probability proportional to the number of data
      # points that already traversed the path
      branch_point <- names(phylobase::descendants(tree_kept_phylo4, root_child, type="child"))
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
          div_t <- div_time(dist_to_root_child, m_v, c, c_order = c_order, theta, alpha)
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
          edge_len <- phylobase::edgeLength(tree_kept_phylo4)[phylobase::getEdge(tree_kept_phylo4, root_child)]
          dist_to_root_child <- start_time + edge_len
        }
      }
    }
  }
}



#' Attach a subtree to a given DDT at a randomly selected location
#' @param subtree subtree to attach to tree_kept
#' @param tree_kept the tree to be attached to
#' @param detach_div_time divergence time of subtree when it was extracted from the original tree
#' @param pa_detach_node_label label of the parent node of the detached node
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#' @param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0. For general multifurcating tree from a Pitman-Yor process,
#'    specify positive values to alpha and theta. It is, however, recommended using alpha = 
#'    theta = 0 in inference because multifurcating trees have not been tested rigorously.
#' @importFrom ape bind.tree
#' @family sample trees
attach_subtree <- function(subtree, tree_kept, detach_div_time, pa_detach_node_label, c,
                           c_order = 1, theta=0.0, alpha=0.0){
  # phylobase::phylobase.options(singleton="ok")
  # plot(subtree, show.node.label = TRUE)
  # if the subtree contains a single node
  if(length(subtree) == 1){
    # print("this1")
    # select a point on the tree to attach subtree to
    new_attach_point <- reattach_point(tree_kept, c, c_order, theta, alpha)
    tree_kept_phylo4 <- phylobase::phylo4(tree_kept, check.node.labels="keep")
    # tree_kept_phylo42 <- phylo4(tree_kept, check.node.labels="drop")
    # if the new attach time is after than the detach time in the original, then resample the attach time until it is before
    while(new_attach_point$div_time >= detach_div_time){
      new_attach_point <- reattach_point(tree_kept, c=c, c_order=c_order, theta=theta, alpha=alpha)
    }
    # attach the subtree to the given tree as a leaf branch
    new_phylo4 <- phylobase::phylo4( add_leaf_branch(tree_old = tree_kept, div_t = new_attach_point$div_time, new_leaf_label = subtree,
                                          where = getNode(tree_kept_phylo4, new_attach_point$root_child),
                                          position = nodeHeight(tree_kept_phylo4, new_attach_point$root_child, "root") - new_attach_point$div_time) )
    # add label to the new attach point
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))] <- pa_detach_node_label
    
    return(list(new_phylo4 = new_phylo4,
                attach_root = new_attach_point$root_node,
                attach_to = new_attach_point$root_child,
                new_div_time = new_attach_point$div_time))
  }else if( phylobase::nTips(tree_kept) == 1 ){ # if the subtree has only one leaf
    # sample a new divergence time
    div_t <- div_time(t_u = 0, m_v = 1, c = c, c_order = c_order, theta = theta, alpha = alpha)
    # if the new attach time is after than the detach time in the original, then resample the attach time until it is before
    while(div_t >= detach_div_time){
      div_t <- div_time(t_u = 0, m_v = 1, c = c, c_order = c_order, theta = theta, alpha = alpha)
    }
    # add a root node to the subtree
    new_subtree <- add_root(tree_old = subtree,
                            root_edge_length = detach_div_time - div_t,
                            root_label = pa_detach_node_label,
                            leaf_label = names(rootNode(phylo4(subtree))))
    # add the subtree to the given tree
    new_phylo4 <- phylo4(ape::bind.tree(tree_kept, new_subtree, where = 1, position = 1 - div_t))
    # add label to the new attach point
    nodeLabels(new_phylo4)[is.na(nodeLabels(new_phylo4))] <- pa_detach_node_label
    return(list(new_phylo4 = new_phylo4,
                attach_root = tree_kept$node.label,
                attach_to = tree_kept$tip.label,
                new_div_time = div_t))
  } else{
    # select a point on the tree to attach subtree to
    new_attach_point <- reattach_point(tree_kept, c=c, c_order=c_order, theta=theta, alpha=alpha)
    # print(new_attach_point)
    tree_kept_phylo4 <- phylo4(tree_kept, check.node.labels="keep")
    while (new_attach_point$div_time >= detach_div_time){
      new_attach_point <- reattach_point(tree_kept, c=c, c_order=c_order, theta=theta, alpha=alpha)
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



#' Calculate proposal likelihood
#' @description Given an old tree, propose a new tree and calculate the original
#'  and proposal tree likelihood in the DDT process
#' @param old_tree_phylo4 the old "phylo4" object
#' @param tree_kept the remaining "phylo" tree after detachment
#' @param old_detach_pa_div_time a number in (0, 1) indicating the divergence time 
#'  of the detached node on the old tree
#' @param old_pa_detach_node_label a character label of the parent of the detached node 
#'  on the old tree
#' @param old_detach_node_label a character label of the detached node 
#'  on the old tree
#' @param new_div_time a number in (0, 1) indicating the divergence time at which
#'  the detached subtree will be re-attached on the proposal tree
#' @param new_attach_root,new_attach_to a character label of the starting and ending
#'  nodes of the branch on the proposal tree, which the detached subtree will be re-attached to
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#' @return a list containing the following elements:
#' \describe{
#' \item{`q_new`}{a "phylo" tree detached from the input tree}
#' \item{`q_old`}{the remaining "phylo" tree after detachment}
#' }
proposal_log_prob <- function(old_tree_phylo4, tree_kept, old_detach_pa_div_time, old_pa_detach_node_label, old_detach_node_label,
                              new_div_time, new_attach_root, new_attach_to, c, c_order=1){
  # phylobase::phylobase.options(singleton="ok")
  if (c_order) {
    A <- a_t_one_cum
    a_t <- a_t_one
  } else{
    A <- a_t_two_cum
    a_t <- a_t_two
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


#'Sample divergence function parameter c for a(t) = c / (1-t) through Gibbs sampler
#'@param shape0 shape of the Inverse-Gamma prior
#'@param rate0 rate of the Inverse-Gamma prior
#'@param tree_structure a data.frame containing the divergence times and number of
#'    data points to the left and right branches of internal nodes on the tree
#'@return a numeric value of the newly sampled c
sample_c_one <- function(shape0, rate0, tree_structure){
  shape <- shape0 + nrow(tree_structure)
  rate <- rate0 - sum(( sapply(tree_structure$m_v-1, function(x) H_n(x)) -
                          sapply(tree_structure$l-1,function(x) H_n(x)) -
                          sapply(tree_structure$r-1,function(x) H_n(x)) ) *
                        (log(1 - tree_structure$div_times)))
  return(rgamma(1, shape = shape, rate = rate))
}

#'Sample divergence function parameter c for a(t) = c / (1-t)^2 through Gibbs sampler
#'@param shape0 shape of the Inverse-Gamma prior
#'@param rate0 rate of the Inverse-Gamma prior
#'@param tree_structure a data.frame containing the divergence times and number of
#'    data points to the left and right branches of internal nodes on the tree
#'@return a numeric value of the newly sampled c
sample_c_two <- function(shape0, rate0, tree_structure){
  shape <- shape0 + nrow(tree_structure)
  rate <- rate0 - sum(( sapply(tree_structure$m_v-1, function(x) H_n(x)) -
                          sapply(tree_structure$l-1,function(x) H_n(x)) -
                          sapply(tree_structure$r-1,function(x) H_n(x)) ) *
                        tree_structure$div_times / (1 - tree_structure$div_times) )
  return(rgamma(1, shape = shape, rate = rate))
}



#' Sample a new tree topology using Metropolis-Hastings through randomly
#'    detaching and re-attaching subtrees
#'@param tree_phylo4d_old a phylo4d object of tree from the previous iteration
#'@param Sigma_by_group a vector of diffusion variances of G groups from the previous iteration
#'@param item_membership_list a vector of G elements, each indicating the number of
#'    items in this group
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#'@param tree_structure_old a data.frame of tree structure from the previous iteration. Each row
#'    contains information of an internal node, including divergence times, number of data points
#'    traveling through the left and right branches
#'@param dist_mat_old a list of leaf covariance matrix from the previous iteration. The list has length
#'    G, the number of item groups
#'@return a numeric vector of G elements, each being the newly sampled variance
#'    of the latent location of this group
sample_tree_topology <- function(tree_phylo4d_old, Sigma_by_group, item_membership_list, c, c_order=1,
                                 tree_structure_old = NULL, dist_mat_old = NULL){
  # number of leaves / classes
  K <- nTips(tree_phylo4d_old)
  old_tree <- extractTree(tree_phylo4d_old)
  # randomly detach a subtree
  random_detach <- random_detach_subtree(tree_phylo4 = old_tree)
  # randomly attach the detached subtree
  proposal_newtree <- suppressWarnings({
    attach_subtree(subtree = random_detach$tree_detached, tree_kept = random_detach$tree_kept,
                   detach_div_time = random_detach$detach_div_time,
                   pa_detach_node_label = random_detach$pa_detach_node_label,
                   c_order = c_order, c = c, theta=0, alpha = 0)
  })
  # while(nTips(proposal_newtree$new_phylo4) != K){
  # need to avoid extremely short branch in the new proposal tree to avoid degeneracy
  # in the covariance matrix
  while(proposal_newtree$new_div_time > 0.95){
    proposal_newtree <- suppressWarnings({attach_subtree(subtree = random_detach$tree_detached, tree_kept = random_detach$tree_kept,
                                       detach_div_time = random_detach$detach_div_time,
                                       pa_detach_node_label = random_detach$pa_detach_node_label,
                                       c_order = c_order, c = c, theta=0, alpha = 0)})
  }
  # compute log likelihood of the old and proposal tree divergence times
  q_prob <- suppressWarnings({proposal_log_prob(old_tree_phylo4 = old_tree, tree_kept = random_detach$tree_kept,
                              old_detach_pa_div_time = random_detach$pa_div_time,
                              old_pa_detach_node_label = random_detach$pa_detach_node_label,
                              old_detach_node_label = random_detach$detach_node_label,
                              new_div_time = proposal_newtree$new_div_time,
                              new_attach_root = proposal_newtree$attach_root,
                              new_attach_to = proposal_newtree$attach_to,
                              c = c, c_order = c_order)})
  # convert to phylo4d object
  new_node_order <- match(proposal_newtree$new_phylo4@label, tree_phylo4d_old@label) # c(paste0("v", 1:K), paste0("u", 1:K))
  # reorder rows of data to from leaf to internal nodes
  #tree_phylo4d_old@data drops all NA rows automatically
  tdata(tree_phylo4d_old)[is.na(tdata(tree_phylo4d_old))] <- 0
  tree_phylo4d_old@data <- tree_phylo4d_old@data[as.character(1:(2*K)),]
  data_new_order <- tree_phylo4d_old@data[as.character(new_node_order),]
  rownames(data_new_order) <- 1:(2*K)
  
  # if (nrow(tree_phylo4d_old@data) == (2*K)){ # all internal node parameters are not NA
  #   tree_phylo4d_old@data <- tree_phylo4d_old@data[as.character(1:(2*K)),]
  #   data_new_order <- tree_phylo4d_old@data[as.character(new_node_order),]
  #   rownames(data_new_order) <- 1:(2*K)
  # } else{ #tree_phylo4d_old@data drops all NA rows automatically
  #   # J <- ncol(tree_phylo4d_old@data)
  #   tdata(tree_phylo4d_old)[is.na(tdata(tree_phylo4d_old))] <- 0
  #   # tr_phylo4d <- phylo4d(tree_phylo4d_old, node.data = matrix(0, nrow = K, ncol = J), merge.data = TRUE)
  #   # tree_phylo4d_old@data <- tree_phylo4d_old@data[as.character(1:(1+K)),]
  #   
  # }
  new_phylo4d <- suppressWarnings({phylo4d(proposal_newtree$new_phylo4, all.data=data_new_order)})

  # compute log likelihood of the old and proposal tree locatioins
  old_tree_info <- logllk_ddt(c = c, c_order, Sigma_by_group = Sigma_by_group, tree_phylo4d = tree_phylo4d_old,
                              item_membership_list = item_membership_list,
                              tree_structure_old = tree_structure_old, dist_mat_old = dist_mat_old)
  new_tree_info <- logllk_ddt(c = c, c_order, Sigma_by_group = Sigma_by_group, tree_phylo4d = new_phylo4d,
                              item_membership_list = item_membership_list)
  
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







