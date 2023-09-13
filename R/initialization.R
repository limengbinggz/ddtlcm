#' Estimate an initial response profile from latent class model using poLCA()
#'@param K number of latent classes
#'@param data a N by J observed binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@param ... optional arguments for the poLCA function
#'@family initialization functions
#'@return a K by J probability matrix, the k,j-th entry being the response probability to item
#'    j of an individual in class k
initialize_poLCA <- function(K, data, ...){
  variable_names <- colnames(data)
  fm <- as.formula(paste("cbind(", paste(variable_names, collapse= ","), ") ~ 1"))
  # data[,] <- as.integer(data[,])
  model <- poLCA(fm, data.frame(data+1), nclass=K, ...)
  # model <- poLCA(fm, data.frame(data+1), nclass=K, maxiter=100,
  #                tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE)
  response_prob_simple_lcm <- matrix(0, nrow = K, ncol = ncol(data))
  for (j in 1:ncol(data)) {
    response_prob_simple_lcm[,j] <- model$probs[[variable_names[j]]][,2]
  }
  colnames(response_prob_simple_lcm) <- variable_names
  return(list(response_probs = response_prob_simple_lcm,
              class_probability = model$P,
              class_assignments = model$predclass,
              model_poLCA = model))
}


#' Provide a random initial response profile based on latent class mode
#'@param K number of latent classes
#'@param data a N by J observed binary matrix, where the i,j-th element is the response
#'    of item j for individual i
#'@return a K by J probability matrix, the k,j-th entry being the response probability to item
#'    j of an individual in class k
initialize_randomLCM <- function(K, data){
  J <- ncol(data)
  N <- nrow(data)
  # randomly sample initial response probabilities
  response_probs <- matrix(runif(K*J), nrow=K, ncol=J)

  class_probability <- rep(1/K, K)
  class_assignments <- sample(1:K, N, prob = class_probability, replace = TRUE)
  return(list(response_probs = response_probs,
              class_probability = class_probability,
              class_assignments = class_assignments))
}




#' Estimate an initial binary tree on latent classes using hclust()
#'@param leaf_data a K by J matrix of logit(theta_{kj})
#'@param c hyparameter of divergence function a(t)
#'@param c_order equals 1 (default) or 2 to choose divergence function
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
#'@importFrom ape as.phylo bind.tree
#'@family initialization functions
#'@return phylo4d object of tree topology
initialize_hclust <- function(leaf_data, c, c_order=1, method_dist = "euclidean", method_hclust = "ward.D",
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
    root_edge_length <- div_time(0, 1, c, c_order, alpha, theta)
    # print(root_edge_length)
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




#' Initialize the MH-within-Gibbs algorithm for DDT-LCM
#'@param K number of classes (integer)
#'@param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#'@param item_membership_list a list of G elements, where the g-th element contains the column
#'  indices of `data` corresponding to items in major group g
#'@param c hyparameter of divergence function a(t)
#'@param c_order equals 1 (default) or 2 to choose divergence function
#'  a(t) = c/(1-t) or c/(1-t)^2.
#'@param method_lcm a character. If `random` (default), the initial LCM parameters will be random values. 
#'  If `poLCA`, the initial LCM parameters will be EM algorithm estimates from the \code{poLCA} function.
#'@param method_dist string specifying the distance measure to be used in dist().
#'    This must be one of "euclidean" (defaults), "maximum", "manhattan", "canberra", "binary" or "minkowski".
#'    Any unambiguous substring can be given.
#'@param method_hclust string specifying the distance measure to be used in hclust().
#'    This should be (an unambiguous abbreviation of) one of "ward.D" (defaults), "ward.D2", "single",
#'    "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or
#'    "centroid" (= UPGMC).
#'@param method_add_root string specifying the method to add the initial branch to the tree
#'    output from hclust(). This should be one of "min_cor" (the absolute value of the minimum
#'    between-class correlation; default) or "sample_ddt" (randomly sample a small divergence time from
#'    the DDT process with c = 100)
#'@param fixed_initials a named list of fixed initial values, including
#' the initial values for tree ("phylo4d"), response_prob, class_probability, class_assignments,
#' Sigma_by_group, and c. Default is NULL. See 
#'@param fixed_priors a named list of fixed prior hyperparameters, including the 
#' the Gamma prior for `c`, inverse-Gamma prior for `sigma_g^2`, and Dirichlet prior
#' for `pi`. Moreover, we allow for a type III generalized logistic distribution such 
#' that f(`eta`; a_pg) = `theta`. This becomes a standard logistic distribution when a_pg = 1. See
#' Dalla Valle, L., Leisen, F., Rossini, L., & Zhu, W. (2021). A Pólya–Gamma sampler for a 
#' generalized logistic regression. Journal of Statistical Computation and Simulation, 91(14), 2899-2916.
#' An example input list is 
#' `list("shape_c" = 1, "rate_c" = 1, "shape_sigma" = rep(2, G), "rate_sigma" = rep(2, G), "a_pg" = 1.0)`, where
#' `G` is the number of major item groups. Default is NULL. 
#'@param alpha,theta hyparameter of branching probability a(t) Gamma(m-alpha) / Gamma(m+1+theta)
#'    For DDT, alpha = theta = 0
#'@param ... optional arguments for the poLCA function
#'@return phylo4d object of tree topology
#'@family initialization functions
#'@seealso [ddtlcm_fit()]
#'@export
#'@examples
#' # load the MAP tree structure obtained from the real HCHS/SOL data
#' data(data_synthetic)
#' # extract elements into the global environment
#' list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv()) 
#' K <- 3
#' G <- length(item_membership_list)
#' fixed_initials <- list("shape_c" = 2, "rate_c" = 2)
#' fixed_priors <- list("rate_sigma" = rep(3, G))
#' initials <- initialize(K, data = response_matrix, item_membership_list, 
#'   c=1, c_order=1, fixed_initials = fixed_initials, fixed_priors = fixed_priors)
initialize <- function(K, data, item_membership_list, c=1, c_order=1,
                       method_lcm = "random", #c("poLCA", "random"),
                       method_dist = "euclidean", method_hclust = "ward.D",
                       method_add_root = "min_cor", 
                       fixed_initials = list(), fixed_priors = list(),
                       alpha=0, theta=0, ...){
  # phylobase::phylobase.options(singleton="ok")
  G <- length(item_membership_list)
  J <- length(unlist(item_membership_list))

  ### step 1: initialize LCM
  if (length(fixed_initials) != 6){
    if (method_lcm == "poLCA"){
      # run classical LCM on the responses
      simple_lcm <- initialize_poLCA(K=K, data, ...)
      # simple_lcm <- initialize_poLCA(K=K, data, maxiter=100,
      #                 tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE)
    } else if (method_lcm == "random"){
      simple_lcm <- initialize_randomLCM(K=K, data)
    } else {
      stop("Initialization method for LCM should be chosen from 'poLCA' or 'random'.")
    }
    response_prob <- simple_lcm$response_probs
    class_probability <- simple_lcm$class_probability
    class_assignments <- simple_lcm$class_assignments
    leaf_data <- logit(response_prob)
  } else {
    initials <- fixed_initials
  }
  if (!is.null(fixed_initials$response_prob)){
    response_prob <- fixed_initials$response_prob
    leaf_data <- logit(response_prob)
  }
  if (!is.null(fixed_initials$class_probability)){
    class_probability <- fixed_initials$class_probability / sum(fixed_initials$class_probability)
  }
  if (!is.null(fixed_initials$class_assignments)){
    class_assignments <- fixed_initials$class_assignments
  }


  ### step 2: estimate c using MLE
  if (!is.null(fixed_initials$tree_phylo4d)){
    tree_phylo4d <- fixed_initials$tree_phylo4d
    init_info <- logllk_ddt(1, c_order, Sigma_by_group = rep(1, G), tree_phylo4d=tree_phylo4d,
                            item_membership_list,
                            tree_structure_old = NULL, dist_mat_old = NULL)
    c_part <- log( 1 - init_info$tree_structure$div_times) *
      exp(init_info$tree_structure$logllk_tree_topology)
    c_hclust <- -(K-1) / sum(c_part)
  } else{
    # get the information of the initial tree by adding a small initial segment
    initial_tree_hclust <- initialize_hclust(leaf_data, 100, c_order, method_dist,
                                             method_hclust, method_add_root, alpha=0, theta=0)
    # add leaf and root data to the tree
    root_location <- rep(0, J)
    # add root to the leaf data
    tree_data <- rbind(leaf_data, root_location)
    # need to add rows corresponding to internal nodes
    tree_data <- rbind(tree_data, matrix(NA, nrow = K-1, ncol = ncol(tree_data)))
    rownames(tree_data) <- c(paste0("v", 1:K), paste0("u", 1:K))
    colnames(tree_data) <- c(paste0("x", 1:ncol(tree_data)))
    initial_tree_hclust <- addData(initial_tree_hclust, all.data = tree_data )
    init_info <- logllk_ddt(1, c_order, Sigma_by_group = rep(1, G), tree_phylo4d=initial_tree_hclust,
                            item_membership_list,
                            tree_structure_old = NULL, dist_mat_old = NULL)
    c_part <- log( 1 - init_info$tree_structure$div_times) *
      exp(init_info$tree_structure$logllk_tree_topology)
    c_hclust <- -(K-1) / sum(c_part)
    
    ### step 3: run hierarchical clustering on logit response profiles
    # options(warn = -1)
    initial_tree_hclust <- initialize_hclust(leaf_data, c_hclust, c_order, method_dist,
                                             method_hclust, method_add_root, alpha=0, theta=0)
    
    # step 4: add leaf and root data to the tree
    # let the root location be the population logit mean
    root_location <- rep(0, J)
    # add root to the leaf data
    tree_data <- rbind(leaf_data, root_location)
    # need to add rows corresponding to internal nodes
    tree_data <- rbind(tree_data, matrix(NA, nrow = K-1, ncol = ncol(tree_data)))
    rownames(tree_data) <- c(paste0("v", 1:K), paste0("u", 1:K))
    colnames(tree_data) <- c(paste0("x", 1:ncol(tree_data)))
    tree_phylo4d <- addData(initial_tree_hclust, all.data = tree_data )
  }
  
  if (!is.null(fixed_initials$c)){
    c_hclust <- fixed_initials$c
  }
  

  ### step 5: estimate Sigma_by_group
  if (!is.null(fixed_initials$Sigma_by_group)){
    Sigma_by_group <- fixed_initials$Sigma_by_group
  } else{
    Sigma_by_group <- rep(0, G)
    for (g in 1:G) {
      Sigma_by_group[g] <- (sd(leaf_data[, item_membership_list[[g]] ]))^2
    }
  }  
  
  ### return initial values
  initials = list(
    tree_phylo4d = tree_phylo4d,
    response_prob = response_prob,
    class_probability = class_probability,
    class_assignments = class_assignments,
    Sigma_by_group = Sigma_by_group,
    c = c_hclust
  )

  priors <- list(
    shape_sigma = rep(2, G),
    rate_sigma = rep(2, G),
    prior_dirichlet = rep(5, K),
    shape_c = 1,
    rate_c = 1,
    a_pg = 1.0
  )
  if (!is.null(fixed_priors$shape_sigma)){
    priors$shape_sigma <- fixed_priors$shape_sigma
  }
  if (!is.null(fixed_priors$rate_sigma)){
    priors$rate_sigma <- fixed_priors$rate_sigma
  }
  if (!is.null(fixed_priors$prior_dirichlet)){
    priors$prior_dirichlet <- fixed_priors$prior_dirichlet
  }
  if (!is.null(fixed_priors$shape_c)){
    priors$shape_c <- fixed_priors$shape_c
  }
  if (!is.null(fixed_priors$rate_c)){
    priors$rate_c <- fixed_priors$rate_c
  }
  if (!is.null(fixed_priors$a_pg)){
    priors$a_pg <- fixed_priors$a_pg
  }
  
  if ("simple_lcm" %in% ls()) {
    return(list(initials = initials, priors = priors, model_lcm = simple_lcm))
  } else{
    return(list(initials = initials, priors = priors))
  }
}



