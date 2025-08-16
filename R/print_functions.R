#' Print out setup of a ddt_lcm model
#' @param x a "ddt_lcm" object
#' @param \dots	Further arguments passed to each method
#' @rdname print.ddt_lcm
#' @method print ddt_lcm
#' @family ddt_lcm results
#' @export
#' @return NULL
#' @examples
#' data(result_diet_1000iters)
#' print(result_diet_1000iters)
print.ddt_lcm <- function(x, ...){
  # total number of items
  J <- length(unlist(x$setting$item_membership_list))
  # number of major item groups
  G <- length(x$setting$item_membership_list)
  N <- nrow(x$data)
  K <- x$setting$K
  total_iters <- length(x$loglikelihood)

  cat("\n---------------------------------------------\n")
  cat("DDT-LCM with K =", K, "latent classes run on", N, "observations and", J, "items in", G, "major groups.\n")
  cat(total_iters, "iterations of posterior samples drawn.")
  cat("\n---------------------------------------------\n")
}


#' Print out summary of a ddt_lcm model
#' @param x a "summary.ddt_lcm" object
#' @param digits integer indicating the number of decimal places (round) to be used.
#' @param \dots Further arguments passed to each method
#' @method print summary.ddt_lcm
#' @family ddt_lcm results
#' @export
#' @return NULL
#' @examples
#' data(result_diet_1000iters)
#' burnin <- 500
#' summarized_result <- summary(result_diet_1000iters, burnin, relabel = TRUE, be_quiet = TRUE)
#' print(summarized_result)
print.summary.ddt_lcm <- function(x, digits = 3L, ...){
  # total number of items
  J <- length(unlist(x$setting$item_membership_list))
  # number of major item groups
  G <- length(x$setting$item_membership_list)
  N <- nrow(x$data)
  K <- x$setting$K
  total_iters <- length(x$loglikelihood)

  cat("Conditional item response (column) probabilities of positive responses, by outcome variable, for each class (row):\n\n")
  item_names <- colnames(x$data)
  response_probs_mean <- matrix(x$response_probs_summary[,"Mean"], nrow = K)
  rownames(response_probs_mean) <- paste0("Class", 1:K)
  colnames(response_probs_mean) <- item_names
  print(round(response_probs_mean, digits))
  cat("\n\n")

  cat("Estimated class probabilities: \n", round(x$class_probs_summary[,"Mean"], digits), "\n\n")

  cat("Maximum log-likelihood of full model:", round(x$max_llk_full, digits), "\n \n")
  cat("Maximum log-likelihood of LCM only:", round(x$max_llk_lcm, digits), "\n \n")
}




