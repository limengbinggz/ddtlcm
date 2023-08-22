#' Generic function
#' @param object The object to summarize
#' @param \dots Further arguments passed to each method
#' @export
summary <- function(object, ...) UseMethod('summary')



#' Summarize the output of a ddt_lcm model
#' @param result a "ddt_lcm" object
#' @param burnin number of samples to discard from the posterior chain as burn-ins. Default is 3000.
#' @param relabel If TRUE, perform post-hoc label switching using the Equivalence Classes 
#'  Representatives (ECR) method to solve non-identifiability issue in mixture models. If FALSE,
#'  no label switching algorithm will be performed. 
#' @param be_quiet If TRUE, do not print information during summarization. If FALSE, print label 
#'  switching information and model summary.
#' @rdname summary.ddt_lcm
#' @method summary ddt_lcm
#' @importFrom label.switching label.switching
#' @family ddt_lcm results
#' @return \code{NULL}
#' @export
summary.ddt_lcm <- function(object, burnin = 3000, relabel = T, be_quiet = F, ...){
  if (!inherits(object, "ddt_lcm")){
    stop("result should be a class ddt_lcm object.")
  }
  total_iters <- length(object$loglikelihood)
  item_membership_list <- object$setting$item_membership_list
  # total number of items
  J <- length(unlist(item_membership_list))
  # number of major item groups
  G <- length(item_membership_list)
  Sigma_by_group <- rep(1, G)
  K <- object$setting$K
  # num_items_per_group <- unlist(lapply(1:G, function(x) rep(x, length(item_membership_list[[x]]))))
  num_items_per_group <- unlist(lapply(item_membership_list, length))
  
  ## get the MAP tree estimate
  map_index <- which.max(object$loglikelihood[-(1:burnin)])
  if (length(object$tree_samples$tree_list) > 1){
    tree_map <- object$tree_samples$tree_list[[map_index+burnin]]
    # get the information of the initial tree
    tree_info <- logllk_ddt(1, object$controls$c_order, Sigma_by_group, tree_map, item_membership_list,
                            tree_structure_old = NULL, dist_mat_old = NULL)
    tree_Sigma <- tree_info$dist_mat$dist_mat_g
  } else{
    tree_map <- NULL
    tree_Sigma <- NULL
  }
  
  
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
  Sigma_summary <- posterior_summary(x=t(object$Sigma_by_group_samples[,(burnin+1):total_iters,drop=F]), var_names)
  # divergence parameter
  c_summary <- posterior_summary(object$c_samples[(burnin+1):total_iters,drop=F], "c")
  
  response_probs_samples <- object$response_probs_samples[(burnin+1):total_iters,,]
  # response_probs_samples <- array(t(object$response_probs_samples[,(burnin+1):total_iters]),
  #                                 dim = c(total_iters-burnin, K, J))
  class_probs_samples <- t(object$class_probs_samples)[(burnin+1):total_iters,,drop=F]
  
  ### if we need to relabel
  if (relabel) {
    map_index <- which.max(object$loglikelihood_lcm[(burnin+1):total_iters]) + burnin
    
    # require(label.switching)
    # switch labels. m = # of iterations
    # sink("NUL")
    quiet(
      ls_lcm <- label.switching(
        method = c("ECR"),
        zpivot = object$Z_samples[,map_index],
        # mxN integer array of the latent allocation vectors generated from an MCMC algorithm
        z = t(object$Z_samples[,(burnin+1):total_iters]),
        K = K,
        # KxJ array containing the parameter that will be used as a pivot
        prapivot = matrix(object$response_probs_samples[map_index,,], nrow=K),
        constraint = 1,
        mcmc = response_probs_samples
      ), be_quiet
    )
    # sink()
    for (iter in 1:(total_iters-burnin)){
      response_probs_samples[iter,,] <- response_probs_samples[iter, ls_lcm$permutations$`ECR`[iter,],]
      class_probs_samples[iter,] <- class_probs_samples[iter, ls_lcm$permutations$`ECR`[iter,]]
    }
  }
  # latent class proportions
  class_probs_summary <- posterior_summary(class_probs_samples, paste0("pi_", 1:K))
  # class response profiles
  ### need to update the extraction of J_g from object ###
  # var_names <- data.frame(k=rep(1:K, each = J), g=rep(1:G, num_items_per_group), j=unlist(sapply(num_items_per_group, seq)))
  var_names <- data.frame(k=rep(1:K, J), g=rep(1:G, num_items_per_group), j = rep(c(unlist(sapply(num_items_per_group, seq))), each = K) )
  var_names <- paste0("theta_", apply(var_names, 1, paste0, collapse = ","))
  response_probs_samples <- matrix(apply(response_probs_samples, 3, c), ncol = K*J)
  response_probs_summary <- posterior_summary(response_probs_samples, var_names)
  max_llk_full <- object$loglikelihood[map_index]
  max_llk_lcm <- object$loglikelihood_lcm[map_index]
  out <- list(tree_map = tree_map, tree_Sigma = tree_Sigma, response_probs_summary = response_probs_summary,
              class_probs_summary = class_probs_summary,
              Sigma_summary = Sigma_summary, c_summary = c_summary,
              max_llk_full = max_llk_full, max_llk_lcm = max_llk_lcm, 
              setting = object$setting, data = object$data)
  class(out) <- "summary.ddt_lcm"
  
  if (!be_quiet){
    print(out)
  }
  return(out)
}
