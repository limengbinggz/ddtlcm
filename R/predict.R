
#' Prediction of class memberships from posterior summaries
#' @description Predict individual class memberships based on posterior summary (point estimates
#'  of model parameters). The predicted class memberships are modal assignments. 
#' @param model_summary a "summary.ddt_lcm" object
#' @param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#' @export
predict.summary.ddt_lcm <- function(model_summary, data){ 
  class_probability <- model_summary$class_probs_summary[,'Mean']
  K <- length(class_probability)
  N <- nrow(data)
  leaf_data <- logit(matrix(model_summary$response_probs_summary[,'Mean'], nrow = K))
  log_probs <- log(class_probability) + leaf_data %*% t(data - 1) + colSums(t(log_expit(leaf_data)))
  # log_probs <- log(class_probability) + leaf_data %*% t(data) - colSums(t(log1p(exp(leaf_data))))
  # substract the max
  log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
  predictive_probs <- apply(log_probs, 2, exp_normalize)
  class_assignments <- rep(0, N)
  # for each individual
  for (i in 1:N) {
    class_assignments[i] <- sample(1:K, 1, prob = predictive_probs[,i]) 
  }
  return(list(class_assignments = class_assignments, predictive_probs = t(predictive_probs)))
}


#' Prediction of class memberships from posterior predictive distributions
#' @description Predict individual class memberships based on posterior predictive distributions.
#'  For each posterior sample, let the class memberships be modal assignments. Then aggregate over
#'  all posterior samples to obtain the most likely assigned classes.
#' @param model a ddt_lcm object
#' @param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#' @param burnin number of samples to discard from the posterior chain as burn-ins. Default is 3000.
#' @export
predict.ddt_lcm <- function(model, data, burnin = 3000){ 
  total_iters <- length(model$loglikelihood)
  K <- model$setting$K
  
  class_assignment_freq <- matrix(0, nrow = N, ncol = K)
  for (iter in (burnin+1):total_iters) {
    class_probability <- model$class_probs_samples[,iter]
    leaf_data <- logit(model$response_probs_samples[iter,,])
    log_probs <- log(class_probability) + leaf_data %*% t(data - 1) + colSums(t(log_expit(leaf_data)))
    log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
    predictive_probs <- apply(log_probs, 2, exp_normalize)
    class_assignments <- apply(predictive_probs, 2, which.max)
    # class_assignments <- rep(0, N)
    # for each individual
    for (i in 1:N) {
      class_assignment_freq[i, class_assignments[i]] <- class_assignment_freq[i, class_assignments[i]] + 1
    }
  }
  class_assignment_freq <- class_assignment_freq / (total_iters-burnin+1)
  class_assignments <- apply(class_assignment_freq, 1, which.max)
  return(list(class_assignments = class_assignments, predictive_probs = class_assignment_freq))
}

