###############################################################
#' MH-within-Gibbs sampler to sample from the full posterior
#' distribution of DDT-LCM
###############################################################

#' @description Use DDT-LCM to estimate latent class and tree on class profiles for 
#'  multivariate binary outcomes.  
#' @param K number of classes (integer)
#' @param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#' @param item_membership_list a list of G elements, where the g-th element contains the column
#'  indices of `data` corresponding to items in major group g, and G is number of major item groups
#' @param total_iters number of posterior samples to collect (integer)
#' 
#' @param initials a named list of initial values of the following parameters:
#' \describe{
#' \item{`tree_phylo4d`}{a phylo4d object. The initial tree have K leaves (labeled as "v1" through "vK"),
#'  1 singleton root node (labeled as "u1"), and K-1 internal nodes (labeled as "u1" through "u{K-1}"). 
#'  The tree also contains parameters for the leaf nodes and the root node (which 
#'  equals 0). The parameters for the internal nodes can be NAs because they will not be used in the
#'  algorithm.}
#' \item{`response_prob`}{a K by J matrix with entries between `0` and `1`. The initial values for the 
#'  item response probabilities. They should equal to the expit-transformed leaf parameters of `tree_phylo4d`.}
#' \item{`class_probability`}{a K-vector with entries between 0 and 1. The initial values for the class 
#'  probabilities. Entries should be nonzero and sum up to 1, or otherwise will be normalized}
#' \item{`class_assignments`}{a N-vector with integer entries from {1, ..., K}. The initial values for
#'  individual class assignments.}
#' \item{`Sigma_by_group`}{a G-vector greater than 0. The initial values for the group-specific diffusion
#'  variances.}
#' \item{`c`}{a value greater than 0. The initial values for the group-specific diffusion
#'  variances.}
#' }
#' Parameters not supplied with initial values will be initialized using the \code{initialize} function
#'  with arguments in `initialize_args`.
#'
#'@param priors a named list of values of hyperparameters of priors. See the function 
#'  \code{initialize} for explanation.
#' \describe{
#' \item{`shape_sigma`}{a G-vector of positive values. The g-th element is the shape parameter for the 
#'  inverse-Gamma prior on diffusion variance parameter sigma_g^2. Default is rep(2, G).}
#' \item{`rate_sigma`}{a G-vector of positive values. Rate parameter. See above. Default is rep(2, G).}
#' \item{`prior_dirichlet`}{a K-vector with entries positive entries. The parameter of the Dirichlet prior
#'  on class probability.}
#' \item{`shape_c`}{a positive value. The shape parameter for the Gamma prior on divergence function 
#'  hyperparameter `c`. Default is 1.}
#' \item{`rate_c`}{a positive value. The rate parameter for `c`. Default is 1.}
#' \item{`a_pg`}{a positive value. The scale parameter for the generalized logistic distribution used in
#'  the augmented Gibbs sampler for leaf parameters. Default is 1, corresponding to the standard logistic
#'  distribution.}
#' }
#'
#'@param controls a named list of control variables. 
#' \describe{
#' \item{`fix_tree`}{a logical. If `TRUE` (default), the tree structure will be sampled in the algorithm. If `FALSE`,
#'  the tree structure will be fixed at the initial input.}
#' \item{`c_order`}{a numeric value. If `1`, the divergence function is a(t) = c/(1-t). If `2`, the divergence 
#' function is a(t) = c/(1-t)^2.}
#' }
#' 
#'@param initialize_args a named list of initialization arguments. See the function 
#'  \code{initialize} for explanation.
#' 
#' @return an object of class "ddt_lcm"; a list containing the following elements:
#' \describe{
#' \item{`tree_samples`}{a list of information of the tree collected from the sampling algorithm, including: 
#'  `accept`: a binary vector where `1` indicates acceptance of the proposal tree and `0` indicates rejection. 
#'  `tree_list`: a list of posterior samples of the tree.
#'  `dist_mat_list`: a list of tree-structured covariance matrices representing the marginal covariances 
#'  among the leaf parameters, integrating out the internal node parameters and all intermediate stochastic paths
#'  in the DDT branching process.}
#' \item{`response_probs_samples`}{a `total_iters` x `K` x `J` array of posterior samples of item response probabilities}
#' \item{`class_probs_samples`}{a `K` x `total_iters` matrix of posterior samples of class probabilities}
#' \item{`Z_samples`}{a `N` x `total_iters` integer matrix of posterior samples of individual class assignments}
#' \item{`Sigma_by_group_samples`}{a `G` x `total_iters` matrix of posterior samples of diffusion variances}
#' \item{`c_samples`}{a `total_iters` vector of posterior samples of divergence function hyperparameter}
#' \item{`loglikelihood`}{a `total_iters` vector of log-likelihoods of the full model}
#' \item{`loglikelihood_lcm`}{a `total_iters` vector of log-likelihoods of the LCM model only}
#' \item{`setting`}{a list of model setup information, including: `K`, `item_membership_list`, and `G`}
#' \item{`controls`}{a list of model controls, including: 
#'  `fix_tree`: FALSE to perform MH sampling of the tree, TRUE to fix the tree at the initial input.
#'  `c_order`: a numeric value of `1` or `2` (see Arguments))}
#' \item{`data`}{the input data matrix}
#' }
#' 
#'@examples
#'# load the MAP tree structure obtained from the real HCHS/SOL data
#'data(data_synthetic)
#'# extract elements into the global environment
#'list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv()) 
#'# run DDT-LCM
#'result <- ddtlcm_fit(K = 3, data = response_matrix, item_membership_list, total_iters = 50)
#'@export
ddtlcm_fit <- function(K, data, item_membership_list, total_iters = 5000, 
                   initials = list(), priors = list(),
                   controls = list(),
                   initialize_args = list(
                     method_lcm = "random", method_dist = "euclidean", 
                     method_hclust = "ward.D", method_add_root = "min_cor",
                     alpha=0, theta=0
                   )
                   ){
  
  ### extract dimensions
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  # number of subjects
  N <- nrow(data)
  if (ncol(data) != J){
    stop("Number of columns in data should be equal to total number of items.")
  }
  item_labels <- paste0("x", 1:J)

  ### initialize parameters
  initialization <- R.utils::doCall(initialize, K = K, data = data, 
                  item_membership_list = item_membership_list, c = 1, c_order=1,
                  fixed_initials = initials, fixed_priors = priors,
                  args = initialize_args)
  initials <- initialization$initials
  priors <- initialization$priors
  
  
  ### extract initial values
  tree_phylo4d <- initials$tree_phylo4d
  if(nTips(tree_phylo4d) != K) stop("The initial tree should have the same number of leaves as the number of latent classes.")
  c <- initials$c
  Sigma_by_group <- initials$Sigma_by_group
  class_assignments <- initials$class_assignments
  class_probability <- initials$class_probability
  tree_phylo4d_old <- tree_phylo4d
  tree_structure_old = NULL
  dist_mat_old = NULL
  leaf_data <- as.matrix(tree_phylo4d@data[names(tipLabels(tree_phylo4d)), item_labels])
  rownames(leaf_data) <- tipLabels(tree_phylo4d)
  leaf_node_labels <- paste0("v", 1:K)
  
  # controls
  if (!is.null(controls$fix_tree)) {
    fix_tree <- controls$fix_tree
  } else{
    fix_tree <- FALSE
  }
  if (!is.null(controls$c_order)) {
    c_order <- controls$c_order
  } else{
    c_order <- 1
  }
  
  if (c_order == 1){
    sample_c <- sample_c_one
  } else if (c_order == 2){
    sample_c <- sample_c_two
  } else {
    stop("c_order must take value 1 or 2.")
  }
  
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

  # initialize storage
  tree_structure_list <- list()
  tree_list <- list()
  dist_mat_list <- list()
  c_samples <- rep(0, total_iters)
  Sigma_by_group_samples <- matrix(NA, length(Sigma_by_group), total_iters)
  # response_probs_samples <- matrix(NA, nrow = K*J, ncol = total_iters)
  response_probs_samples <- array(0, dim = c(total_iters, K, J))
  class_probs_samples <- matrix(0, nrow = K, ncol = total_iters)
  Z_samples <- matrix(0, nrow = N, ncol = total_iters)
  
  # vector of acceptance
  accept <- rep(0, total_iters)
  # model log likelihoods
  logllk_model <- rep(0, total_iters)
  logllk_model_lcm <- rep(0, total_iters)
  
  # get the information of the initial tree
  init_info <- logllk_ddt(c, c_order, Sigma_by_group, tree_phylo4d, item_membership_list,
                          tree_structure_old = NULL, dist_mat_old = NULL)
  # tree_list[[1]] <- tree_phylo4d
  # dist_mat_list[[1]] <- init_info$dist_mat
  leaf_nodes <- tipLabels(tree_phylo4d)
  # logllk_model[1] <- init_info$logllk +
  #   logllk_lcm(data, tree_phylo4d@data[names(leaf_nodes), item_labels],
  #              class_probability, prior_dirichlet, ClassItem, Class_count)
  internal_nodes <- nodeLabels(tree_phylo4d)
  if (fix_tree){
    dist_mat_old <- init_info$dist_mat
    tree_structure_old <- init_info$tree_structure
  }
  old_leaf_labels <- tipLabels(tree_phylo4d)
  
  # set.seed(2022)
  cat("\n## Start posterior sampling ##")
  for(iter in 1:total_iters){
    if (iter %% 500 == 0){
      cat("\n## iteration ", iter, "completed.")
    }
    
    
    ### Sample tree topology
    if (!fix_tree){
      sampled_tree <- sample_tree_topology(tree_phylo4d_old, Sigma_by_group, item_membership_list, c = c, c_order,
                                           tree_structure_old = tree_structure_old, dist_mat_old = dist_mat_old)#NULL
      accept[iter] <- sampled_tree$accept
      logllk_model[iter] <- sampled_tree$logllk_model
      tree_structure_old <- sampled_tree$tree_structure
      dist_mat_old <- sampled_tree$dist_mat
      # to save memory, we only update the tree list if the proposal is accepted; otherwise
      # create a pointer to the old object
      if (sampled_tree$accept | iter == 1) {
        tree_list[[iter]] <- sampled_tree$tree_phylo4d
        dist_mat_list[[iter]] <- sampled_tree$dist_mat
        tree_phylo4d_old <- sampled_tree$tree_phylo4d
        # old_leaf_labels <- tipLabels(tree_phylo4d_old)
      } else{
        tree_list[[iter]] <- tree_list[[iter-1]]
        dist_mat_list[[iter]] <- dist_mat_list[[iter-1]]
      }
    }

    
    ### Sample divergence time hyperparameter and group variance
    c <- sample_c(shape0 = shape_c, rate0 = rate_c, tree_structure = tree_structure_old)
    c_samples[iter] <- c
    Sigma_by_group <- sample_sigmasq(shape0 = shape_sigma, rate0 = rate_sigma, dist_mat = dist_mat_old,
                                     item_membership_list, leaf_data)
    Sigma_by_group_samples[,iter] <- Sigma_by_group
    
    
    ### Sample leaf locations and polya-gamma auxiliary variables
    # if (iter %% 100 == 0){
    sampled <- sample_leaf_locations_pg(item_membership_list, dist_mat_old,
                                        Sigma_by_group, pg_mat, a_pg, auxiliary_mat, 
                                        auxiliary_mat_range, class_assignments)
    # replace old leaf locations with the newly sampled ones
    leaf_data <- sampled$leaf_data
    rownames(leaf_data) <- leaf_node_labels
    colnames(leaf_data) <- item_labels
    new_node_order <- match(tipLabels(tree_phylo4d_old), leaf_node_labels)
    tree_phylo4d_old@data[as.character(1:K) , item_labels] <- leaf_data[new_node_order,]
    pg_mat <- sampled$pg_data
    auxiliary_mat <- sampled$auxiliary_mat

        
    ### Sample individual class assignments Z_i
    class_assignments <- sample_class_assignment(data, leaf_data, a_pg, auxiliary_mat, class_probability)
    Class_count <- tabulate(class_assignments, nbins = K) * 1.0
    ClassItem <- matrix(0, nrow = K, ncol = J)
    for (k in 1:K) {
      class_k <- data[class_assignments == k,,drop = FALSE]
      if (length(class_k) > 0){
        ClassItem[k,] <- colSums(data[class_assignments == k,,drop = FALSE])
      } else{
        ClassItem[k,] <- 0
      }
    }

    ### sample class probability
    class_probability <- c(rdirichlet(1, prior_dirichlet + Class_count))
    response_probs_samples[iter,,] <- expit(leaf_data)
    class_probs_samples[,iter] <- class_probability
    Z_samples[, iter] <- class_assignments
    
    ### compute model loglikelihood
    logllk_model_lcm[iter] <- logllk_lcm(data, leaf_data, class_probability, prior_dirichlet, ClassItem, Class_count)
    logllk_model[iter] <- logllk_model[iter] + logllk_model_lcm[iter]
  }
  cat("\n## Finish posterior sampling\n")
  
  # combine all results in a list
  tree_samples <- list()
  tree_samples[["accept"]] <- accept
  tree_samples[["tree_list"]] <- tree_list
  tree_samples[["dist_mat_list"]] <- dist_mat_list
  setting <- list(K = K, item_membership_list = item_membership_list, G = G)
  controls <- list(fix_tree = fix_tree, c_order = c_order)
  result <- list(
    tree_samples = tree_samples,
    response_probs_samples = response_probs_samples,
    class_probs_samples = class_probs_samples,
    Z_samples = Z_samples,
    Sigma_by_group_samples = Sigma_by_group_samples,
    c_samples = c_samples,
    loglikelihood = logllk_model,
    loglikelihood_lcm = logllk_model_lcm,
    setting = setting,
    controls = controls,
    data = data
  )
  class(result) <- "ddt_lcm"
  return(result)
}







