library(poLCA)

#'@description estimate an initial response profile from latent class model using poLCA()
#'@param K number of latent classes
#'@param response_matrix a N by J observed binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param ... optional arguments for the poLCA function
#'@return a K by J probability matrix, the k,j-th entry being the response probability to item
#'    j of an individual in class k
initialize_poLCA <- function(K, response_matrix, ...){
  variable_names <- colnames(response_matrix)
  fm <- as.formula(paste("cbind(", paste(variable_names, collapse= ","), ") ~ 1"))
  # response_matrix[,] <- as.integer(response_matrix[,])
  model <- poLCA(fm, data.frame(response_matrix+1), nclass=K, ...)
  # model <- poLCA(fm, data.frame(response_matrix+1), nclass=K, maxiter=100,
  #                tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE)
  response_prob_simple_lcm <- matrix(0, nrow = K, ncol = ncol(response_matrix))
  for (j in 1:ncol(response_matrix)) {
    response_prob_simple_lcm[,j] <- model$probs[[variable_names[j]]][,2]
  }
  colnames(response_prob_simple_lcm) <- variable_names
  return(list(response_probs = response_prob_simple_lcm,
              class_probability = model$P,
              class_assignments = model$predclass,
              model_poLCA = model))
}


#'@description provide a random initial response profile based on latent class mode
#'@param K number of latent classes
#'@param response_matrix a N by J observed binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param ... optional arguments for the poLCA function
#'@return a K by J probability matrix, the k,j-th entry being the response probability to item
#'    j of an individual in class k
initialize_randomLCM <- function(K, response_matrix, ...){
  J <- ncol(response_matrix)
  N <- nrow(response_matrix)
  # randomly sample initial response probabilities
  response_probs <- matrix(runif(K*J), nrow=K, ncol=J)

  # leaf_data <- logit(response_probs)
  # # KxN log probabilities: the i-th column is the unnormalized log probabilities of K classes
  # # for indiviudal i
  # class_probability_prior <- rep(1/K, K)
  # log_probs <- log(class_probability_prior) + leaf_data %*% t(response_matrix - 1) + colSums(t(log_expit(leaf_data)))
  # # substract the max
  # log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
  # class_assignments <- rep(0, N)
  # # for each individual
  # for (i in 1:N) {
  #   class_assignments[i] <- sample(1:K, 1, prob = exp_normalize(log_probs[,i]))
  # }
  # return(list(response_probs = response_probs,
  #             class_probability = as.vector(unname(table(class_assignments) / N)),
  #             class_assignments = class_assignments))

  class_probability <- rep(1/K, K)
  class_assignments <- sample(1:K, N, prob = class_probability, replace = T)
  # print(class_probability)
  return(list(response_probs = response_probs,
              class_probability = class_probability,
              class_assignments = class_assignments))
}




#'@description estimate an initial binary tree on latent classes using hclust()
#'@param leaf_data a K by J matrix of logit(theta_{kj})
#'@param c hyparameter of divergence function  a(t) = c / (1-t)
#'@param method_dist string specifying the distance measure to be used in dist().
#'    This must be one of "euclidean", "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'    Any unambiguous substring can be given.
#'@param method_hclust string specifying the distance measure to be used in hclust().
#'    This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'    "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#'    "centroid" (= UPGMC).
#'@param method_add_root string specifying the method to add the initial branch to the tree
#'    output from hclust(). This should be one of "min_cor" (the absolute value of the minimum
#'    between-class correlation) or "sample_ddt" (randomly sample a small divergence time from
#'    the DDT process with a large c = 100)
#'@param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0
#'@param ... optional arguments for the poLCA function
#'@return phylo4d object of tree topology
initialize_hclust <- function(leaf_data, c, is_linear=T, method_dist = "euclidean", method_hclust = "ward.D",
                              method_add_root = "min_cor",
                              alpha=0, theta=0, ...){
  K <- nrow(leaf_data)
  # create distance matrix between classes
  leaf_dist <- dist(leaf_data, method = method_dist)
  # hierarchical clustering
  hclust_phylo <- as.phylo(hclust(leaf_dist, method = method_hclust))
  total_length <- nodeHeight(phylo4d(hclust_phylo), "1", "root")
  # since the resulting tree is not rooted, we need to sample a divergence time
  # for the root node and add to the tree
  if (method_add_root == "min_cor"){
    # use the absolute value of the minimum between-class correlation
    root_edge_length <- abs(min(cor(t(leaf_data))))
  } else if (method_add_root == "sample_ddt"){
    # randomly sample a very small number
    root_edge_length <- div_time(0, 1, c, is_linear, alpha, theta)
    print(root_edge_length)
  }
  # since the total depth rom hclust is not 1, we need to normalize. The tree
  # from hclust will have total depth 1 - root_edge_length
  # cat("edge len: ", hclust_phylo$edge.length)
  hclust_phylo$edge.length <- hclust_phylo$edge.length / total_length * (1 - root_edge_length)
  # cat("\n\nedge len: ", hclust_phylo$edge.length)
  # add root node
  hclust_phylo <- add_root(hclust_phylo, root_edge_length = root_edge_length,
           root_label = "u1", leaf_label = "u2")
  hclust_phylo <- phylo4d(hclust_phylo)
  # add node labels
  tipLabels(hclust_phylo) <- paste0("v", 1:K)
  nodeLabels(hclust_phylo) <- paste0("u", 1:K)
  # nodeHeight(hclust_phylo, "v1", "root")
  # plot(hclust_phylo, show.node=T)
  return(hclust_phylo)
}




initialize <- function(K, response_matrix, num_items_per_group, c=1, is_linear=T,
                       method_lcm = "random", #c("poLCA", ),
                       method_dist = "euclidean", method_hclust = "single",
                       method_add_root = "min_cor",
                       alpha=0, theta=0, ...){
  G <- length(num_items_per_group)

  # step 1: initialize LCM
  # method_lcm = match.arg(method_lcm)
  if (method_lcm == "poLCA"){
    # step 1: run classical LCM on the responses
    simple_lcm <- initialize_poLCA(K=K, response_matrix, ...)
    # simple_lcm <- initialize_poLCA(K=K, response_matrix, maxiter=100,
    #                 tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE)
  } else if (method_lcm == "random"){
    simple_lcm <- initialize_randomLCM(K=K, response_matrix, ...)
    # print(simple_lcm$class_probability)
  } else {
    stop("Initialization method for LCM should be chosen from 'poLCA' or 'random'.")
  }
  response_prob_simple_lcm <- simple_lcm$response_probs
  class_probability_lcm <- simple_lcm$class_probability
  class_assignments_lcm <- simple_lcm$class_assignments
  leaf_data <- logit(response_prob_simple_lcm)

  # step 2: estimate c using MLE
  # plot(initial_tree_hclust, show.node=T, plot.data = F)
  # nodeHeight(initial_tree_hclust, "v1", "root")
  # get the information of the initial tree by adding a small initial segment
  initial_tree_hclust <- initialize_hclust(leaf_data, 100, is_linear, method_dist,
                                           method_hclust, method_add_root, alpha=alpha, theta=theta)
  # add leaf and root data to the tree
  # let the root location be the population logit mean
  # root_location <- logit(colMeans(response_matrix))
  root_location <- rep(0, sum(G))
  # add root to the leaf data
  tree_data <- rbind(leaf_data, root_location)
  # need to add rows corresponding to internal nodes
  tree_data <- rbind(tree_data, matrix(NA, nrow = K-1, ncol = ncol(tree_data)))
  rownames(tree_data) <- c(paste0("v", 1:K), paste0("u", 1:K))
  colnames(tree_data) <- c(paste0("x", 1:ncol(tree_data)))
  initial_tree_hclust <- addData(initial_tree_hclust, all.data = tree_data )
  init_info <- logllk_ddt(1, is_linear, Sigma_by_group = rep(1, G), tree_phylo4d=initial_tree_hclust,
                          num_items_per_group,
                          tree_structure_old = NULL, dist_mat_old = NULL)
  c_part <- log( 1 - init_info$tree_structure$div_times) *
    exp(init_info$tree_structure$logllk_tree_topology)
  c_hclust <- -(K-1) / sum(c_part)

  # step 3: run hierarchical clustering on logit response profiles
  initial_tree_hclust <- initialize_hclust(leaf_data, c_hclust, is_linear, method_dist,
                                           method_hclust, method_add_root, alpha=alpha, theta=theta)

  # step 4: add leaf and root data to the tree
  # let the root location be the population logit mean
  # root_location <- logit(colMeans(response_matrix))
  root_location <- rep(0, sum(G))
  # add root to the leaf data
  tree_data <- rbind(leaf_data, root_location)
  # need to add rows corresponding to internal nodes
  tree_data <- rbind(tree_data, matrix(NA, nrow = K-1, ncol = ncol(tree_data)))
  rownames(tree_data) <- c(paste0("v", 1:K), paste0("u", 1:K))
  colnames(tree_data) <- c(paste0("x", 1:ncol(tree_data)))
  initial_tree_hclust <- addData(initial_tree_hclust, all.data = tree_data )

  # step 5: estimate Sigma_by_group
  # cumulative index of item groups
  cum_J_index <- c(0, cumsum(num_items_per_group))
  Sigma_by_group_lcm <- rep(0, G)
  for (g in 1:G) {
    # Sigma_tree <- init_info$dist_mat$dist_mat_g
    # ratio <- abs(cov(t(leaf_data[,(cum_J_index[g]+1):(cum_J_index[g+1])])) / Sigma_tree)
    # Sigma_by_group_lcm[g] <- mean(ratio[upper.tri(ratio)])
    Sigma_by_group_lcm[g] <- (sd(leaf_data[,(cum_J_index[g]+1):(cum_J_index[g+1])]))^2
  }

  initials = list(
    tree_phylo4d = initial_tree_hclust,
    c = c_hclust,
    Sigma_by_group = Sigma_by_group_lcm,
    class_probability = class_probability_lcm,
    class_assignments = class_assignments_lcm
  )

  priors <- list(
    shape_c = 1,
    rate_c = 1,
    shape_sigma = rep(2, G),
    rate_sigma = rep(2, G),
    # prior_class_probability = prior_class_probability,
    prior_dirichlet = rep(5, K),
    a_pg = 1.0
  )

  return(list(initials = initials, priors = priors, model_lcm = simple_lcm))
}



