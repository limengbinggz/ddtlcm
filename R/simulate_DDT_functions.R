###############################################################
#' Functions to simulate trees and node parameters from a DDT 
#' process.
###############################################################

#' Add a branch to an existing tree according to the branching process of DDT
#' @param tree_old a "phylo" object. The tree (K leaves) to which a new branch will be added.
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#'  a(t) = c/(1-t) or c/(1-t)^2.
#' @param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    Allowable range: 0 <= beta <= 1, and alpha >= -2 beta
#'    For DDT, alpha = theta = 0. For general multifurcating tree from a Pitman-Yor process,
#'    specify positive values to alpha and theta. It is, however, recommended using alpha = 
#'    theta = 0 in inference because multifurcating trees have not been tested rigorously.
#' @importFrom ape Ntip
#' @return a "phylo" object. A tree with K+1 leaves.
add_one_sample <- function(tree_old, c, c_order, theta, alpha){
  tree <- tree_old
  # get information from the existing tree
  # number of leaves from the old tree
  n <- ape::Ntip(tree_old)
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
  root_node <- phylobase::rootNode(tr_phylo4)
  # get the index of the child of the root node
  root_child <- phylobase::descendants(tr_phylo4, root_node, "child")
  # cat("root_child = ", root_child)
  # get the branch length between the root and its child
  root_child_branch_len <- edgeLength(tr_phylo4) [getEdge(tr_phylo4, names(root_child)) ]
  dist_to_root_child <- start_time + root_child_branch_len
  while(TRUE){
    # sample a new divergence time
    div_t <- div_time(t_u = start_time, m_v = cur_n, c = c, c_order = c_order, theta = theta, alpha = alpha)
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
          div_t <- div_time(dist_to_root_child, m_v, c, c_order = c_order, theta, alpha)
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



#' Simulate a tree from a DDT process. Only the tree topology and branch lengths
#'  are simulated, without node parameters.
#' @param K number of leaves (classes) on the tree
#' @param c hyparameter of divergence function a(t)
#' @param c_order equals 1 (default) or 2 to choose divergence function
#'  a(t) = c/(1-t) or c/(1-t)^2.
#' @param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0. For general multifurcating tree from a Pitman-Yor process,
#'    specify positive values to alpha and theta. It is, however, recommended using alpha = 
#'    theta = 0 in inference because multifurcating trees have not been tested rigorously.
#' @importFrom ape read.tree
#' @importFrom Rdpack reprompt
#' @return A class "phylo" tree with K leaves. The leaf nodes are labeled "v1", ..., "vK", 
#'  root node "u1", and internal nodes "u2", ..., "uK". Note that this tree does not contain
#'  any node parameters.
#' @family simulate DDT-LCM data
#' @references {
#'  Knowles, D. A., & Ghahramani, Z. (2014). Pitman yor diffusion trees for bayesian hierarchical 
#'  clustering. IEEE transactions on pattern analysis and machine intelligence, 37(2), 271-289.
#' }
#' @export
#' @examples
#' K <- 6
#' c <- 5
#' c_order <- 1
#' tree1 <- simulate_DDT_tree(K, c, c_order)
#' tree2 <- simulate_DDT_tree(K, c, c_order, alpha = 0.4, theta = 0.1)
#' tree3 <- simulate_DDT_tree(K, c, c_order, alpha = 0.8, theta = 0.1)
simulate_DDT_tree <- function(K, c, c_order = 1, alpha = 0, theta = 0){
  # the first data point does not diverge, and reaches t = 1 directly.
  # the second data point diverges at some time t_2
  div_t <- div_time(t_u = 0.0, m_v = 1.0, c = c, c_order = c_order, theta = theta, alpha = alpha)
  # create the initial tree with two branches
  tree_txt <- paste("((1:", 1-div_t, ",2:",1-div_t,"):", div_t, ");", sep='')
  tree <- read.tree(text = tree_txt)
  # plot(as(tree,"phylo4"))
  
  # Now the third data point and so on ------------------------------------
  for (n in 3:K) {
    tree <- add_one_sample(tree, c, c_order, theta, alpha)
  }
  tree$node.label <- paste0("u", 1:Nnode(tree))
  return(tree)
}


#' Simulate node parameters along a given tree.
#'@param tree_phylo a "phylo" object containing the tree topology and branch lengths
#'@param Sigma_by_group a G-vector greater than 0. The initial values for the group-specific diffusion
#'  variances
#'@param item_membership_list a list of G elements, where the g-th element contains the indices of 
#'  items in major group g
#'@param root_node_location the coordinate of the root node parameter. By default, the node parameter
#'  initiates at the origin so takes value 0. If a value, then the value will be repeated into a 
#'  length J vector. If a vector, it must be of length J.
#'@return A class "phylo4d" tree with K leaves with node parameters. The leaf nodes are labeled "v1", ..., "vK", 
#'  root node "u1", and internal nodes "u2", ..., "uK". 
#'@importFrom ape Ntip Nnode
#'@family simulate DDT-LCM data
#'@export
#'@examples
#' library(ape)
#' tr_txt <- "(((v1:0.25, v2:0.25):0.65, v3:0.9):0.1);"
#' tree <- read.tree(text = tr_txt)
#' tree$node.label <- paste0("u", 1:Nnode(tree))
#' plot(tree, show.node.label = TRUE)
#' # create a list of item membership indices of 7 major groups
#' item_membership_list <- list()
#' num_items_per_group <- c(rep(10, 5), 15, 15)
#' G <- length(num_items_per_group)
#' j <- 0
#' for (g in 1:G) {
#'   item_membership_list[[g]] <- (j+1):(j+num_items_per_group[g])
#'   j <- j+num_items_per_group[g]
#' }
#' # variance of logit response probabilities of items in each group
#' Sigma_by_group <- c(rep(0.6**2, 5), rep(2**2, 2)) #rep(1**2, G)
#' set.seed(50)
#' tree_with_parameter <- simulate_parameter_on_tree(tree, Sigma_by_group, item_membership_list)
simulate_parameter_on_tree <- function(tree_phylo, Sigma_by_group, item_membership_list, 
                                       root_node_location = 0){
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  num_items_per_group <- unlist(lapply(item_membership_list, length))
  if (length(Sigma_by_group) != G){
    stop("Sigma_by_group must be a length G vector, where G is the number of major item groups.")
  }
  if (any(Sigma_by_group <= 0)){
    stop("Sigma_by_group must contain positive values.")
  }
  
  # add labels to all nodes in the tree
  tree_phylo$tip.label <- paste("v", 1:ape::Ntip(tree_phylo), sep="")
  tree_phylo$node.label <- paste("u", 1:ape::Nnode(tree_phylo), sep="")
  tr_phylo4 <- as(tree_phylo, "phylo4")
  # initialize data storage
  location_data_list <- list() 
  # now we start with the root node
  start_node <- rootNode(tr_phylo4)
  if (length(root_node_location) == 1){
    location_data_list[[names(start_node)]] <- rep(root_node_location, J)
  } else if (length(root_node_location) != J){
    stop("The dimension of root location of the DDT should be the same as the number of items J.")
  } else {
    location_data_list[[names(start_node)]] <- root_node_location
  }
  
  while( length(start_node)!=0 ){
    next_from_nonTip <- numeric()
    for (cur_from_idx in seq_along(start_node)){
      from_node <- start_node[cur_from_idx]
      to_nodes <- descendants(tr_phylo4, node=names(from_node), type="children")
      for (i in seq_along(to_nodes)) {
        to_node = to_nodes[i]
        edge_len <- edgeLength(tr_phylo4)[getEdge(tr_phylo4, node=names(to_node), type="descendant")]
        
        Sigma <- rep(Sigma_by_group, num_items_per_group) * edge_len
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
  tr_phylo4d <- phylo4d(tr_phylo4, all.data = location_data, metadata = list(item_membership_list = item_membership_list))
  return(tr_phylo4d)
}




#' Simulate multivariate binary responses from a latent class model
#'@description Generate multivariate binary responses from the following process:
#'    For individual i = 1, ..., N, draw \eqn{Z_i} from Categorical distribution with prior class probability (length K).
#'    For item j = 1, ..., J, given \eqn{Z_i = k}, draw \eqn{Y_{ij}} from Binomial with class-item probability
#'@param N number of individuals
#'@param response_prob a K by J matrix, where the k,j-th element is the response
#'    probability of item j for individuals in class k
#'@param class_probability a length K vector, where the k-th element is the
#'    probability of assigning an individual to class k. It does not have to sum up to 1
#'@return a named list of the following elements:
#'\describe{
#' \item{`response_matrix`}{a K by J matrix with entries between `0` and `1` for the item 
#'  response probabilities.}
#' \item{`class_probability`}{a K-vector with entries between 0 and 1 for the class 
#'  probabilities. Entries should be nonzero and sum up to 1, or otherwise will be normalized}
#' }
#' @family simulate DDT-LCM data
#' @export
#' @examples
#' # number of latent classes
#' K <- 6
#' # number of items
#' J <- 78
#' response_prob <- matrix(runif(K*J), nrow = K)
#' class_probability <- rep(1/K, K)
#' # number of individuals
#' N <- 100
#' response_matrix <- simulate_lcm_response(N, response_prob, class_probability)
simulate_lcm_response <- function(N, response_prob, class_probability){
  # number of classes
  K <- nrow(response_prob)
  # number of items
  J <- ncol(response_prob)
  # initialize binary reponse matrix
  response_matrix <- matrix(0, nrow = N, ncol = J)
  class_assignment <- rep(0, N)
  for (i in 1:N) {
    # class assignment
    Z <- sample(1:K, 1, prob = class_probability)
    class_assignment[i] <- Z
    response_matrix[i,] <- rbinom(J, size = 1, prob = response_prob[Z,])
  }
  colnames(response_matrix) <- paste0("x", 1:J)
  return(list(response_matrix = response_matrix, class_assignment = class_assignment))
}





#' Simulate multivariate binary responses from a latent class model given a tree
#'@description Generate multivariate binary responses from the following process:
#'    For individual i = 1, ..., N,
#'        \eqn{Z_i ~ Categorical_K(prior_class_probability)}
#'        For item j = 1, ..., J,
#'            \eqn{Y_{ij} | Z_i = k ~ Binomial(class_item_probability_{kj})}
#'@param tree_phylo a "phylo" tree with K leaves
#'@param N number of individuals
#'@param class_probability a length K vector, where the k-th element is the
#'  probability of assigning an individual to class k. It does not have to sum up to 1
#'@param item_membership_list a list of G elements, where the g-th element contains the 
#'  column indices of the observed data matrix corresponding to items in major group g
#'@param Sigma_by_group a length-G vector for the posterior mean group-specific diffusion variances.
#'@param root_node_location the coordinate of the root node parameter. By default, the node parameter
#'  initiates at the origin so takes value 0. If a value, then the value will be repeated into a 
#'  length J vector. If a vector, it must be of length J.
#'@param seed_parameter an integer random seed to generate parameters given the tree
#'@param seed_response an integer random seed to generate multivariate binary observations from LCM
#'@return a named list of the following elements:
#'\describe{
#' \item{`tree_with_parameter`}{a "phylo4d" tree with K leaves.}
#' \item{`response_prob`}{a K by J matrix, where the k,j-th element is the response
#'  probability of item j for individuals in class k}
#' \item{`response_matrix`}{a K by J matrix with entries between 0 and 1 for the item 
#'  response probabilities.}
#' \item{`class_probability`}{a K-vector with entries between 0 and 1 for the class 
#'  probabilities. Entries should be nonzero and sum up to 1, or otherwise will be normalized}
#' \item{`class_assignments`}{a N-vector with integer entries from 1, ..., K. The initial values for
#'  individual class assignments.}
#' \item{`Sigma_by_group`}{a G-vector greater than 0. The initial values for the group-specific diffusion
#'  variances.}
#' \item{`c`}{a value greater than 0. The initial values for the group-specific diffusion
#'  variances.}
#' \item{`item_membership_list`}{same as input}
#' }
#' @family simulate DDT-LCM data
#' @export
#' @examples
#' # load the MAP tree structure obtained from the real HCHS/SOL data
#' data(parameter_diet)
#' # unlist the elements into variables in the global environment
#' list2env(setNames(parameter_diet, names(parameter_diet)), envir = globalenv()) 
#' # number of individuals
#' N <- 496
#' # set random seed to generate node parameters given the tree
#' seed_parameter = 1
#' # set random seed to generate multivariate binary observations from LCM
#' seed_response = 1
#' # simulate data given the parameters
#' sim_data <- simulate_lcm_given_tree(tree_phylo, N, 
#'                                     class_probability, item_membership_list, Sigma_by_group, 
#'                                     root_node_location = 0, seed_parameter = 1, seed_response = 1)
simulate_lcm_given_tree <- function(tree_phylo, N, 
            class_probability = 1, item_membership_list, Sigma_by_group = NULL, 
            root_node_location = 0, seed_parameter = 1, seed_response = 1) {

  if (!inherits(tree_phylo, "phylo")){
    stop("Argument `tree_phylo` should be a 'phylo' object.")
  }
  set.seed(seed_parameter)
  K <- phylobase::nTips(tree_phylo)
  # print(K)
  # simulate node parameters along the tree
  tree_phylo4d <- simulate_parameter_on_tree(tree_phylo, Sigma_by_group, item_membership_list, root_node_location)
  
  # simulate multivariate binary observations
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  # extract logit of response profiles by class
  leaf_data <- as.matrix(tree_phylo4d@data[as.character(1:K), 1:J])
  response_probs <- expit(leaf_data)
  set.seed(seed_response)
  # simulate multivariate responses
  if (length(class_probability) == 1){
    class_probability <- rep(class_probability, K)
  }
  sim_responses <- simulate_lcm_response(N, response_probs, class_probability)
  # N x J binary response matrix
  response_matrix <- sim_responses$response_matrix
  # true class assignment
  class_assignment_true <- sim_responses$class_assignment
  
  # save simulated data as a list
  sim_data <- list()
  sim_data[["tree_with_parameter"]] <- tree_phylo4d
  sim_data[["class_probability"]] <- class_probability
  sim_data[["response_prob"]] <- response_probs
  sim_data[["response_matrix"]] <- response_matrix
  sim_data[["class_assignment"]] <- class_assignment_true
  sim_data[["Sigma_by_group"]] <- Sigma_by_group
  sim_data[["item_membership_list"]] <- item_membership_list
 
  return(sim_data) 
}


