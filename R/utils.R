
# divergence function a(t) = c / (1-t)
a_t_one <- function(c, t) c/(1-t)
# cumulative hazard function for a(t) = c / (1-t)
A_t_one <- function(c, t) -c * log(1-t)
#' inverse divergence function for a(t) = c / (1-t)
A_t_inv_one <- function(c, y) 1.0 - exp(- y/c)

# divergence function a(t) = c / (1-t)^2
a_t_two <- function(c, t) c/(1-t)^2
# cumulative hazard function for a(t) = c / (1-t)^2
A_t_two <- function(c, t) -c + c / (1.0-t)
#' inverse divergence function for a(t) = c / (1-t)^2
A_t_inv_two <- function(c, y) 1.0 - c / (c+y)

# sample divergence time on an edge uv previously traversed by m(v) data points
div_time <- function(t_u, m_v, c, c_order = 1, alpha = 0, theta = 0){
  u <- runif(1)
  if (c_order == 1){
    x = A_t_one(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u)
    return (A_t_inv_one(c, x))
  } else if (c_order == 2){
    x = A_t_two(c, t_u) - exp( lgamma(m_v+1.0+theta) - lgamma(m_v-alpha) ) * log(1-u);
    return (A_t_inv_two(c, x))
  } else {
    stop("c_order must take value 1 or 2.")
  }
}


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
                        # where = as.numeric(substr(names(where), 2, nchar(names(where)))),
                        where = where,
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
  tree_new <- bind.tree(tree_old, branch, 
                        where = where
                        # where = as.numeric(substr(names(where), 2, nchar(names(where))))
                        )
  return(tree_new)
}



add_root <- function(tree_old, root_edge_length, root_label, leaf_label){
  leaf_branch <- list(edge = matrix(c(2,1), nrow = 1),
                      node.label = root_label,
                      tip.label = leaf_label,
                      edge.length = root_edge_length,
                      Nnode = 1)
  class(leaf_branch) <- "phylo"
  tree_new <- bind.tree(leaf_branch, tree_old, where=1)
  return(tree_new)
}


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


#' The expit function 
#' @description The expit function: f(x) = exp(x) / (1+exp(x)), computed
#'  in a way to avoid numerical underflow.
#' @param x a value or a numeric vector between 0 and 1 (exclusive)
#' @export
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


#' The logistic function
#' @description The logit function: f(x) = log(x / (1/x)). Large absolute values of 
#'  x will be truncated to +/- 5 after logit transformation according to its sign.
#' @param x a value or a numeric vector between 0 and 1 (exclusive)
#' @export
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


#' suppress print from cat()
quiet <- function(x, be_quiet=T) { 
  if (be_quiet){
    sink(tempfile()) 
    on.exit(sink()) 
    invisible(force(x)) 
  }else{
    x
  }
} 
