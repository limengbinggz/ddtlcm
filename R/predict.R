#' Prediction of class memberships from posterior summaries
#' @description Predict individual class memberships based on posterior summary (point estimates
#'  of model parameters). The predicted class memberships are modal assignments.
#' @param object a "summary.ddt_lcm" object
#' @param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#' @param \dots Further arguments passed to each method
#' @export
#' @return a list of the following named elements:
#' \describe{
#' \item{`class_assignments`}{an integer vector of individual predicted class memberships
#'  taking values in 1, ..., K}
#' \item{`predictive_probs`}{a N x K matrix of probabilities, where the (i,k)-th element
#'  is the probability that the i-th individual is predicted to belong to class k.}
#' }
#' @examples
#' data(result_diet_1000iters)
#' burnin <- 500
#' summarized_result <- summary(result_diet_1000iters, burnin, relabel = TRUE, be_quiet = TRUE)
#' predicted <- predict(summarized_result, result_diet_1000iters$data)
predict.summary.ddt_lcm <- function(object, data, ...){
  class_probability <- object$class_probs_summary[,'Mean']
  K <- length(class_probability)
  N <- nrow(data)
  leaf_data <- logit(matrix(object$response_probs_summary[,'Mean'], nrow = K))
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
#' @param object a ddt_lcm object
#' @param data an NxJ matrix of multivariate binary responses, where
#'   N is the number of individuals, and J is the number of granular items
#' @param burnin number of samples to discard from the posterior chain as burn-ins. Default is 3000.
#' @param \dots Further arguments passed to each method
#' @export
#' @return a list of the following named elements:
#' \describe{
#' \item{`class_assignments`}{an integer vector of individual predicted class memberships
#'  taking values in 1, ..., K}
#' \item{`predictive_probs`}{a N x K matrix of probabilities, where the (i,k)-th element
#'  is the probability that the i-th individual is predicted to belong to class k.}
#' }
#' @examples
#' data(result_diet_1000iters)
#' burnin <- 500
#' predicted <- predict(result_diet_1000iters, result_diet_1000iters$data, burnin)
predict.ddt_lcm <- function(object, data, burnin = 3000, ...){
  total_iters <- length(object$loglikelihood)
  K <- object$setting$K
  N <- nrow(data)

  class_assignment_freq <- matrix(0, nrow = N, ncol = K)
  for (iter in (burnin+1):total_iters) {
    class_probability <- object$class_probs_samples[,iter]
    leaf_data <- logit(object$response_probs_samples[iter,,])
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

