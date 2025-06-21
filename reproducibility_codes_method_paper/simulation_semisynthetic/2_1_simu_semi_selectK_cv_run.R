#######################################################
#' Simulate data along tree (K=6) estimated from real
#' data. Generate different sets of theta along the given tree.
#' For each set of theta, generate multiple multivariate
#' binary responses.
#' Select the best number of latent classes using predictive 
#' likelihood from 5-fold cross-validation.
#######################################################
library(ggtree)
library(clue)
library(aricode)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

source("functions/utils.R")
source("functions/simulate_DDT_functions.R")
source("functions/loglikehoods.R")
source("functions/MH_GLogit.R")
source("functions/initialization.R")
source("functions/summary_functions.R")
options(warn = -1)

# args = commandArgs(trailingOnly=TRUE)
# N <- as.integer(args[1])
# seed_theta = as.integer(args[2])
# seed_Y = as.integer(args[3])
# K_candidate = as.integer(args[4])
# fold_num = as.integer(args[5]) # index of the training fold

## specify parameters -----------------------------
N <- 400
seed_theta <- 18
seed_Y <- 5
K_candidate <- 6
fold_num <- 1

# the number of classes K
K_true <- 6
load("hchs_simulated_data/tree_real_nodata.RData")
num_items_per_group <- J_g
G <- length(J_g)
item_group_membership <- rep(1:G, J_g)
# total number of items
J <- sum(J_g) #78
class_probability_true <- class_probability
Sigma_by_group <- c(2.3, 1, 1, 2.3, 2.3, 1, 2.3)^2


### Simulate data -----------------------------------------------------------
# location of root node in the DDT; can use sample average
root_node_location0 <- rep(0, J)

# we use divergence function a(t) = c / (1-t)
set.seed(seed_theta)
tree_phylo <- as.phylo(extractTree(tree))
tree_with_profiles <- gen_location_with_selector(tree_phylo, Sigma_by_group, item_group_membership,
                                                 root_node_location = root_node_location0)

# extract logit of response profiles by class
leaf_data <- as.matrix(tree_with_profiles@data[as.character(1:K_true), 1:J])
class_item_probability <- expit(leaf_data)
# simulate multivariate responses
set.seed(seed_Y)
sim_responses <- gen_response(N, class_item_probability, class_probability)
# N x J binary response matrix
response_matrix <- sim_responses$response_matrix
# true class assignment
class_assignment_true <- sim_responses$class_assignment

# save simulated data as a list
sim_data <- list()
sim_data[["tree_with_data"]] <- tree_with_selector
sim_data[["class_probability"]] <- class_probability
sim_data[["response_prob"]] <- class_item_probability
sim_data[["response_matrix"]] <- response_matrix
sim_data[["class_assignment"]] <- class_assignment_true
sim_data[["Sigma_by_group"]] <- Sigma_by_group



# split individuals into 5 folds
seed <- fold_num + 9246*seed_Y + 109*seed_theta
set.seed(seed)
training_index <- sample(1:N, as.integer(N/5*4), replace = FALSE)
training_data <- response_matrix[training_index,]


# Estimate the initial tree from LCM --------------------------------------
c0 <- 5
c_init <- rgamma(1, c0, 1)
# initialize randomly
initialization <- initialize(K_candidate, training_data, num_items_per_group = J_g, c = c_init,
                             method_lcm = "random", method_dist = "euclidean", #"poLCA"
                             method_hclust = "single", method_add_root = "min_cor",
                             alpha=0, theta=0, maxiter=100,
                             tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE)

# # alternatively, we may initialize using poLCA
# # poLCA may fail. If this is the case, use another seed
# noError = F
# initialization_num <- 1
# repeat{
#   if (noError == T) break
#   set.seed(initialization_num)
#   initialization <- try(initialize(K, training_data, num_items_per_group = J_g, c = c_init,
#                                    method_lcm = "poLCA", method_dist = "euclidean", #"poLCA"
#                                    method_hclust = "single", method_add_root = "min_cor",
#                                    alpha=0, theta=0, maxiter=100,
#                                    tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE), silent=TRUE)
#   if ('try-error' %in% class(initialization)) initialization_num <- sample.int(9000, size=1)
#   else noError = T
# }

priors <- initialization$priors
initials <- initialization$initials

total_iters <- 8000
burnin <- 5000
controls <- list(
  fix_tree = FALSE,
  is_linear = TRUE
)
system.time({
  result <- ddt_lcm_mcmc(K = K_candidate, data = training_data, num_items_per_group = J_g,
                         initials = initials, priors = priors, 
                         total_iters = total_iters, controls = controls)
})


### calculate the predictive likelihood
testing_data <- response_matrix[-training_index,]

## posterior summary
s <- summary(result, burnin = burnin, relabel = T) #8000
# get class prevalences
class_probs <- s$class_probs_summary[,"Mean"]
## get class profiles 
leaf_data <- logit(matrix(s$response_probs_summary[,'Mean'], nrow = K_candidate))

# predict individual class memberships
log_probs <- log(class_probs) + leaf_data %*% t(testing_data - 1) + colSums(t(log_expit(leaf_data)))
# substract the max
log_probs <- t(t(log_probs) - apply(log_probs, 2, max))
predictive_probs <- apply(log_probs, 2, exp_normalize)
class_prediction <- apply(predictive_probs, 2, which.max)

### calculate predictive likelihood following the LCM posteriro summary ---------------
predictive_llk <- 0
for (k in 1:K_candidate) {
  predictive_llk = predictive_llk + sum(testing_data[class_prediction == k,,drop=F] %*% leaf_data[k,]) 
}
predictive_llk = predictive_llk + 
  tabulate(class_prediction, K_candidate) %*% colSums(t(-log1p(exp(leaf_data)))) +
  sum(tabulate(class_prediction, K_candidate) * log(class_probs))
cat("predictive llk based on summary:", predictive_llk)


### compute predictive likelihood from the posterior distribution ---------------
calculate_predict_llk <- function(iter){
  log_class_probs <- log(result$class_probs_samples[,iter])
  response_probs <- t(matrix(result$response_probs_samples[,iter], nrow = K_candidate))
  # loglikelihood of LCM
  logllk_lcm_sample <- 0
  for (i in 1:nrow(testing_data)) {
    logllk_lcm_sample <- logllk_lcm_sample +
      logSumExp(log_class_probs + colSums(testing_data[i,] * log(response_probs) + (1 - testing_data[i,]) * log(1 - response_probs)))
  }
  return(logllk_lcm_sample)
}
predict_llk_dist <- mclapply((burnin + 1):total_iters, calculate_predict_llk, mc.cores = 2L)
predict_llk_dist <- unlist(predict_llk_dist)
predictive_llk_dist_mean <- mean(predict_llk_dist)
cat("mean predictive llk based on distribution:", predictive_llk_dist_mean)


### save data ------------------------------------------
save_dir <- "simulation_semisynthetic/data/2_semi_selectK"
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/N", N, "_J", J)
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/seedtheta", seed_theta, "_seedY", seed_Y)
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
file_name <- paste0(save_dir, "/K", K_candidate, "_fold", fold_num, ".RData")
save(seed_theta, seed_Y, result, K_true,
     K_candidate, fold_num, predictive_llk, training_index, 
     predictive_llk, predict_llk_dist, predictive_llk_dist_mean, file = file_name)
