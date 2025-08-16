

library(phylobase)
library(ape)
library(data.table)
# library(ggplot2)
library(matrixStats)


# divergence function a(t) = c / (1-t)
a_t_linear <- function(c, t) c/(1-t)
# cumulative hazard function for a(t) = c / (1-t)
A_t_linear <- function(c, t) -c * log(1-t)
#' inverse divergence function for a(t) = c / (1-t)
A_t_inv_linear <- function(c, y) 1.0 - exp(- y/c)

# divergence function a(t) = c / (1-t)^2
a_t_quadratic <- function(c, t) c/(1-t)^2
# cumulative hazard function for a(t) = c / (1-t)^2
A_t_quadratic <- function(c, t) -c + c / (1.0-t)
#' inverse divergence function for a(t) = c / (1-t)^2
A_t_inv_quadratic <- function(c, y) 1.0 - c / (c+y)

# sample divergence time on an edge uv previously traversed by m(v) data points
div_time <- function(t_u, m_v, c, is_linear, alpha=0, theta=0){
  u <- runif(1)
  if (is_linear){
    x = A_t_linear(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u)
    return (A_t_inv_linear(c, x))
  } else {
    x = A_t_quadratic(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u);
    return (A_t_inv_quadratic(c, x));
  }
}

#' a_t <- function(c, t) c/(1-t)
#'
#' #' cumulative hazard function for a(t) = c / (1-t)
#' A <- function(c, t) -c * log(1-t)
#'
#' #' inverse divergence function for a(t) = c / (1-t)
#' A_inv <- function(c, y) 1.0 - exp(- y/c)
#'
#' #' sample divergence time on an edge uv previously traversed by m(v) data points
#' div_time <- function(t_u, m_v, c, alpha, theta){
#'   u <- runif(1)
#'   x <- A(c, t_u) - exp( lgamma(m_v+1+theta) - lgamma(m_v-alpha) ) * log(1-u)
#'   A_inv(c, x)
#' }

#' #' compute trace of a matrix
#' tr <- function(mat) sum(diag(mat))

#' #' check whether a matrix is invertible
#' check_invertible <- function(mat){
#'   out <- qr(mat)
#'   cond <- out$rank == nrow(mat)
#'   return(cond)
#' }

# compute normalized probabilities: exp(x_i) / sum_j exp(x_j)
exp_normalize <- function(x){
  return(exp(x - logSumExp(x)))
}


#' sample multivariate normal using precision matrix
#' from x ~ N(Q^{-1}a, Q^{-1}), where Q^{-1} is the precision matrix
draw_mnorm <- function(precision_mat, precision_a_vec){
  U <- chol(precision_mat)
  b <- rnorm(nrow(precision_mat))
  backsolve(U, backsolve(U, precision_a_vec, transpose = TRUE) + b)
}


#' expit function: f(x) = exp(x) / (1+exp(x))
expit <- function(x){
  out <- x
  positives <- which(x >= 0)
  if (length(positives) > 0){
    # print("1,")
    out[positives] <- 1 / (1 + exp(-x[positives]))
    # return(1 / (1 + exp(-x)))
  }
  if (length(positives) < length(x)) {
    # print("2,")
    negatives <- which(x < 0)
    out[negatives] <- exp(x[negatives]) / (1 + exp(x[negatives]))
    # return(exp(x) / (1 + exp(x)))
  }
  return(out)
}


logit <- function(x){
  y <- log(x / (1-x))
  large <- which(abs(y) > 5)
  if (length(large) > 0){
    y[large] <- sign(y[large]) * 5
  }
  return(y)
}


log_expit <- function(x){
  out <- x
  idx <- which(x < -33.3)
  if (length(idx) > 0){
    out[idx] <- x[idx]
  }
  idx <- which(x >= -33 & x < -18)
  if (length(idx) > 0){
    out[idx] <- x[idx] - exp(x[idx])
  }
  idx <- which((x >= -18) & (x < 37))
  if (length(idx) > 0){
    out[idx] <- -log1p(exp(-x[idx]))
  }
  return(out)
}
