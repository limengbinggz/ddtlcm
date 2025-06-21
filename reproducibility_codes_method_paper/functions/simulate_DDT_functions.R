
library(phylobase)
library(ape)
library(data.table)
library(ggplot2)

# Simulate DDT ------------------------------------------------------------

#' #' cumulative hazard function for a(t) = c / (1-t)
#' A <- function(c, t) -c * log(1-t)
#'
#' #' inverse divergence function for a(t) = c / (1-t)
#' A_inv <- function(y) 1.0 - exp(- y/c)
#'
#' #' sample divergence time on an edge uv previously traversed by m(v) data points
#' div_time <- function(t_u, m_v, c, alpha, theta){
#'   u <- runif(1)
#'   x <- A(c, t_u) - exp( lgamma(m_v+1+theta) - lgamma(m_v-alpha) ) * log(1-u)
#'   A_inv(x)
#' }

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
                        where = as.numeric(substr(names(where), 2, nchar(names(where)))),
                        position = position)
  return(tree_new)
}


#' add a leaf branch to an existing tree tree_old to make a multichotomus branch
#' div_t: divergence time of the new branch
#' new_leaf_label: the label of the newly added leaf
#' where: node name of to which node in the existing tree the new leaf branch should connect to
add_multichotomous_tip <- function(tree_old, div_t, new_leaf_label, where){
  branch <- list(edge = matrix(c(2,1), nrow = 1),
                 tip.label = new_leaf_label,
                 edge.length = 1 - div_t,
                 Nnode = 1)
  class(branch) <- "phylo"
  tree_new <- bind.tree(tree_old, branch, where = as.numeric(substr(names(where), 2, nchar(names(where)))))
  return(tree_new)
}



#' simulate DDT
sim_DDT_data <- function(Nsamples, c, is_linear, alpha, theta){
  #' the first data point does not diverge, and reaches t = 1 directly.
  #' the second data point diverges at some time t_2
  div_t <- div_time(t_u = 0.0, m_v = 1.0, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
  # create the initial tree with two branches
  tree_txt <- paste("((1:", 1-div_t, ",2:",1-div_t,"):", div_t, ");", sep='')
  tree <- read.tree(text = tree_txt)
  # plot(as(tree,"phylo4"))

  #' Now the third data point and so on ------------------------------------
  for (n in 3:Nsamples) {
    tree <- add_one_sample(tree, c, is_linear, theta, alpha)
  }
  tree$node.label <- paste0("u", 1:Nnode(tree))
  return(tree)
}


add_one_sample <- function(tree_old, c, is_linear, theta, alpha){
  tree <- tree_old
  # get information from the existing tree
  # number of leaves from the old tree
  n <- Ntip(tree_old)
  # add labels to the leaf nodes
  tree$tip.label <- paste("v", 1:n, sep="")
  # add labels to the internal nodes
  tree$node.label <- paste("v", (n+1):(n + Nnode(tree)), sep="")
  # convert to phylo4 object
  tr_phylo4 <- phylo4(tree, check.node.labels = "keep")
  # current number of leaves
  cur_n <- n
  # root node starts from time 0
  start_time = 0
  # get the root node index
  root_node <- rootNode(tr_phylo4)
  # get the index of the child of the root node
  root_child <- phylobase::descendants(tr_phylo4, root_node, "child")
  # cat("root_child = ", root_child)
  # get the branch length between the root and its child
  root_child_branch_len <- edgeLength(tr_phylo4) [getEdge(tr_phylo4, names(root_child)) ]
  dist_to_root_child <- start_time + root_child_branch_len
  while(TRUE){
    # sample a new divergence time
    div_t <- div_time(t_u = start_time, m_v = cur_n, c = c, is_linear = is_linear, theta = theta, alpha = alpha)
    # if t2 < t1, then diverge at t2
    if (div_t < dist_to_root_child){
      new_leaf_label <- paste("v", n+1, sep="")
      tree <- add_leaf_branch(tree, div_t, new_leaf_label,
                            where = root_child,
                            position = dist_to_root_child - div_t)
      # plot(test)
      return(tree)
    }else{
      #' if t2 > t1, then select which path to take, with probability proportional to the number of data
      #' points that already traversed the path
      branch_point <- phylobase::descendants(tr_phylo4, root_child, type="child")
      # print(paste("branch_point: ", names(branch_point)))
      # get the number of branches from this branch point
      K <- length(branch_point)
      # get the number of samples which previously took each branch
      n_k <- unlist(lapply(phylobase::descendants(tr_phylo4, names(branch_point), type="tip"), length))
      # compute the selection probabilities
      probs <- c(n_k - alpha, theta + alpha*K) / (sum(n_k) + theta)
      # select which branch to follow, or to create a new branch
      path_idx_sel <- which.max(rmultinom(1, 1, probs))
      if (path_idx_sel == length(probs)){
        # if we select to create a new branch, the resulting tree will be multichotomous at this branch point
        new_leaf_label <- paste("v", n+1, sep="")
        # print(paste("\nroot_child: ", root_child, "\n"))
        tree <- add_multichotomous_tip(tree, dist_to_root_child, new_leaf_label,
                                     where = root_child)
        return(tree)

      } else{
        # if we choose one of the existing branchs, then we need to sample the divergence time again
        selected_node <- branch_point[path_idx_sel]
        # get the number of data points which previously travel the path
        m_v <- n_k[path_idx_sel]
        if( m_v==1 ){
          # sample the divergence time again
          div_t <- div_time(dist_to_root_child, m_v, c, is_linear = is_linear, theta, alpha)
          # if only one data point has traversed along this path, then we may simply diverge at the new time
          new_leaf_label <- paste("v", n+1, sep="")
          tree <- add_leaf_branch(tree, div_t, new_leaf_label,
                                where = selected_node,
                                position = 1 - div_t)
          # plot(test)
          return(tree)

        } else{
          # if more than one data point has  traversed along this path, then we need to repeat the above
          # process again: starting from the branch point, follow the sampled path, until reaching a path that
          # has been traveled by only one data point
          start_time <- dist_to_root_child
          cur_n <- m_v
          root_node <- root_child
          root_child <- selected_node
          edge_len <- edgeLength(tr_phylo4)[getEdge(tr_phylo4, names(root_child))]
          dist_to_root_child <- start_time + edge_len

        }
      }
    }
  }

}



#' tree: a phylo object for the tree structure
#' sigma_sq: variance of the data in Brownian motion
#' dimension: the dimension of each data point
gen_location_general <- function(tree, sigma_sq, dimension){
  # add labels to all nodes in the tree
  tree$tip.label <- paste("v", 1:Ntip(tree), sep="")
  # tree$node.label <- paste("u", ( Ntip(tree)+1 ):( Ntip(tree)+Nnode(tree) ), sep="")
  tree$node.label <- paste("u", 1:Nnode(tree), sep="")
  tr_phylo4 <- as(tree, "phylo4")
  # initialize data storage
  location_data_list <- list() #vector("list", length = Nnode(tree) + Ntip(tree))
  # now we start with the root node
  start_node <- rootNode(tr_phylo4)
  location_data_list[[names(start_node)]] <- rep(0, dimension)
  # names(location_data_list[[names(start_node)]]) <- names(start_node)
  # colnames(location_data) <- paste0("x", 1:dimension)
  Num_tip_located <- 0
  # base_chol_lt <- chol_PD(sigma_sq)
  # sigma_sq<-base_chol_lt$final_mat
  # base_chol<-base_chol_lt$chol_mat

  chol_prec <- chol(solve(sigma_sq))
  #' generate multivariate normal distribution with precision
  #' mu: mean vector
  #' chol_prec: cholesky decomposition of the precision matrix
  rmultinon_precision <- function(mu, chol_prec){
    b = rnorm(nrow(chol_prec))
    backsolve(chol_prec, backsolve(chol_prec, mu, transpose = TRUE) + b)
  }

  while( length(start_node)!=0 ){
    # cat("start_node =", start_node, "\n")
    next_from_nonTip <- numeric()
    for (cur_from_idx in 1:length(start_node)){
      from_node <- start_node[cur_from_idx]
      # cat("from_node =", from_node, "\n")
      to_nodes <- descendants(tr_phylo4, node=names(from_node), type="children")
      for (i in 1:length(to_nodes)) {
        to_node = to_nodes[i]
        edge_len <- edgeLength(tr_phylo4)[getEdge(tr_phylo4, node=names(to_node), type="descendant")]
        # from_location <- as.matrix(location_data[names(from_node),], nrow=length(from_node))
        to_location <- rmultinon_precision(mu = location_data_list[[names(from_node)]], chol_prec / sqrt(edge_len))
        location_data_list[[names(to_node)]] <- to_location
        next_from_nonTip <- c(next_from_nonTip, to_node[!(names(to_node) %in% tipLabels(tr_phylo4))])
        # cat("next_from_nonTip =", next_from_nonTip, "\n")
      }
    }
    start_node <- next_from_nonTip
  }

  location_data <- do.call(rbind, location_data_list)
  colnames(location_data) <- paste0("x", 1:dimension)
  tr_phylo4d <- phylo4d(tr_phylo4, all.data = location_data)
  tr_phylo4d <- addData(tr_phylo4d,
                        tip.data = t(matrix(rep(sigma_sq[upper.tri(sigma_sq, diag=T)], Ntip(tree)),
                                            nrow = ((1+dimension)*dimension/2))) )
  return(tr_phylo4d)
}



#' tree: a phylo object for the tree structure
#' Sigma_by_group: a vector of length G, where Sigma_g is the variance of the data of items in group g in Brownian motion
#' dimension: the dimension of each data point
gen_location_with_selector <- function(tree, Sigma_by_group, item_group_membership, 
                                       root_node_location = NULL){
  # get the number of groups
  G <- length(unique(item_group_membership))
  # number of items in total
  J <- length(item_group_membership)
  J_g <- c(table(item_group_membership))

  # add labels to all nodes in the tree
  tree$tip.label <- paste("v", 1:Ntip(tree), sep="")
  tree$node.label <- paste("u", 1:Nnode(tree), sep="")
  tr_phylo4 <- as(tree, "phylo4")
  # initialize data storage
  location_data_list <- list() #vector("list", length = Nnode(tree) + Ntip(tree))
  # now we start with the root node
  start_node <- rootNode(tr_phylo4)
  if (is.null(root_node_location)){
    location_data_list[[names(start_node)]] <- rep(0, J)
  } else if (length(root_node_location) != J){
    stop("The dimension of root location of the DDT should be the same as the number of items J.")
  } else {
    location_data_list[[names(start_node)]] <- root_node_location
  }

  while( length(start_node)!=0 ){
    next_from_nonTip <- numeric()
    for (cur_from_idx in 1:length(start_node)){
      from_node <- start_node[cur_from_idx]
      to_nodes <- descendants(tr_phylo4, node=names(from_node), type="children")
      for (i in 1:length(to_nodes)) {
        to_node = to_nodes[i]
        edge_len <- edgeLength(tr_phylo4)[getEdge(tr_phylo4, node=names(to_node), type="descendant")]

        Sigma <- rep(Sigma_by_group, J_g) * edge_len
        # cat("location_data_list[[names(from_node)]] =", location_data_list[[names(from_node)]], "\n")
        to_location <- location_data_list[[names(from_node)]]
        to_location <- rnorm(n = length(Sigma), mean = location_data_list[[names(from_node)]], sd = sqrt(Sigma))

        # cat("to_node =", names(to_node), ", with location =", to_location, "\n\n")
        location_data_list[[names(to_node)]] <- to_location
        next_from_nonTip <- c(next_from_nonTip, to_node[!(names(to_node) %in% tipLabels(tr_phylo4))])
      }
    }
    start_node <- next_from_nonTip
  }

  location_data <- do.call(rbind, location_data_list)
  colnames(location_data) <- paste0("x", 1:J)
  tr_phylo4d <- phylo4d(tr_phylo4, all.data = location_data, metadata = list(item_group_membership = item_group_membership))
  return(tr_phylo4d)
}




#'@description Simulate multivariate responses from a latent class model
#'    For individual i = 1, ..., N,
#'        Z_i ~ Categorical_K(prior_class_probability)
#'        For item j = 1, ..., J,
#'            Y_{ij} | Z_i = k ~ Binomial(class_item_probability_{kj})
#'@param N number of individuals
#'@param class_item_probability a K by J matrix, where the k,j-th element is the response
#'    probability of item j for individuals in class k
#'@param prior_class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@return a N by J binary matrix, where the i,j-th element is the response
#'    of item j for individual i
gen_response <- function(N, class_item_probability, prior_class_probability){
  # number of classes
  K <- nrow(class_item_probability)
  # number of items
  J <- ncol(class_item_probability)
  # initialize binary reponse matrix
  response_matrix <- matrix(0, nrow = N, ncol = J)
  class_assignment <- rep(0, N)
  for (i in 1:N) {
    # class assignment
    Z <- sample(1:K, 1, prob = prior_class_probability)
    class_assignment[i] <- Z
    response_matrix[i,] <- rbinom(J, size = 1, prob = class_item_probability[Z,])
  }
  colnames(response_matrix) <- paste0("x", 1:J)
  return(list(response_matrix = response_matrix, class_assignment = class_assignment))
}





#' plot the tree structure and node locations from a DDT, without knowing the granular diffusion
#' process
#' tree: a phylo object containing the tree structure
#' tree_with_locations: a phylo4 object containing the data point locations
plot_DDT_skeleton <- function(tree, tree_with_locations){
  require(ggraph)
  require(igraph)
  require(adephylo)
  tree_igraph <- as.igraph(tree)
  match(tree_with_locations@label, names(V(tree_igraph)))
  V(tree_igraph)$class <- names(V(tree_igraph))
  p <- ggraph(tree_igraph, layout = "dendrogram") +
    geom_edge_link() +
    geom_node_point() +
    geom_node_label(aes(label=class))
  p_positions <- layer_data(p, i = 2L)
  p_data <- p$data
  # get distance from node to root
  root_node <- rootNode(phylo4d(tree))
  dist_to_root <- dist.nodes(tree)[,root_node]
  match_idx <- match(p_data$name, tree_with_locations@label)
  p$data$x <- dist_to_root[match_idx]
  p$data$y <- tree_with_locations@data[match_idx,]$x1
  p

}


plot_response_profiles_by_node <- function(tree, nodes_to_plot, nrow = 2){
  # find the node index in the tree
  node_index <- match(nodes_to_plot, tree@label)
  dat <- tree@data[as.character(node_index),]
  # get the number of items
  J <- length(grep("x", colnames(tree@data)))
  item_group_membership <- tree@metadata$item_group_membership
  node_labels <- nodeLabels(tree)

  plot_data <- vector("list", length = length(nodes_to_plot))
  for (nn in 1:length(nodes_to_plot)){
    node <- nodes_to_plot[nn]
    plot_data_node <- data.table(node_name = node,
                                 eta = unlist(dat[nn, 1:J]),
                                 group_membership = item_group_membership,
                                 item = 1:J)
    if (node %in% node_labels){# if an internal node
      plot_data_node$group_selector <-  paste0(node, ": (", paste0(unlist(dat[nn, (J+1):ncol(dat)]), collapse = ", "), ")")
    } else{
      plot_data_node$group_selector <-  paste0(node, ": NA")
    }
    plot_data_node$theta <- expit(plot_data_node$eta)
    plot_data[[nn]] <- plot_data_node
  }
  plot_data <- do.call(rbind, plot_data)
  plot_data$item <- factor(plot_data$item, levels = unique(plot_data$item))
  plot_data$group_membership <- factor(plot_data$group_membership, levels = unique(plot_data$group_membership))
  plot_data$group_selector <- factor(plot_data$group_selector, levels = unique(plot_data$group_selector))
  # plot_data$node_name <- factor(plot_data$node_name, levels = nodes_to_plot)

  ggplot(plot_data, aes(item, theta)) +
    geom_col(aes(fill = group_membership), position = "identity") +
    facet_wrap(~group_selector, nrow = nrow) +
    labs(x = "Item", y = "Response probability", fill = "Group") +
    theme_bw()
}


