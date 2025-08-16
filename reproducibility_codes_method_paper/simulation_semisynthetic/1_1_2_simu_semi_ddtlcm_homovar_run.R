#######################################################
#' Simulate data along tree (K=6) estimated from real
#' data. Generate different sets of theta along the given tree.
#' For each set of theta, generate multiple multivariate
#' binary responses.
#' Independent normal priors on the logistic-transformed
#' item responsibilities. Homogeneous variance parameters.
#' Z_i ~ Cat(\pi)
#' P(Y_ij = 1 \mid Z_i = k) = \theta_{kj}
#' logistic(\theta_{kj}) ~ N(0, \sigma^2), for any item j 
#######################################################

library(ggtree)
library(clue)
library(aricode)
library(label.switching)
# curr_dir <- "reproducibility_codes"
# setwd(curr_dir)

source("functions/utils.R")
source("functions/simulate_DDT_functions.R")
source("functions/loglikehoods.R")
source("functions/MH_GLogit.R")
source("functions/initialization.R")
source("functions/summary_functions.R")
options(warn = -1)
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

N <- as.integer(args[1])
seed_theta = as.integer(args[2])
seed_Y = as.integer(args[3])

# ## specify parameters -----------------------------
# seed_theta <- 1
# seed_Y <- 2
# N=400


K <- 6
is_linear <- "TRUE"
load("hchs_simulated_data/tree_real_nodata.RData")

num_items_per_group <- J_g
G <- length(J_g)
item_group_membership <- rep(1:G, J_g)
# total number of items
J <- sum(J_g) #78
class_probability_true <- class_probability
# group-specific diffusion variances
Sigma_by_group <- c(2.3, 1, 1, 2.3, 2.3, 1, 2.3)^2

### Simulate data -----------------------------------------------------------
# location of root node in the DDT; can use sample average
root_node_location0 <- rep(0, J)

is_linear <- ifelse(is_linear == "TRUE", T, F)
set.seed(seed_theta)
tree_phylo <- as.phylo(extractTree(tree))
tree_with_profiles <- gen_location_with_selector(tree_phylo, Sigma_by_group, item_group_membership,
                                                 root_node_location = root_node_location0)

# extract logit of response profiles by class
leaf_data <- as.matrix(tree_with_profiles@data[as.character(1:K), 1:J])
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
sim_data[["tree_with_data"]] <- tree_with_profiles
sim_data[["class_probability"]] <- class_probability
sim_data[["response_prob"]] <- class_item_probability
sim_data[["response_matrix"]] <- response_matrix
sim_data[["class_assignment"]] <- class_assignment_true
sim_data[["Sigma_by_group"]] <- Sigma_by_group


# Estimate the initial tree from LCM --------------------------------------
c_init <- rgamma(1, 5, 1)
# initialize randomly
# note that num_items_per_group = J instead of J_g
initialization <- initialize(K, response_matrix, num_items_per_group = J, c = c_init,
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
#   initialization <- try(initialize(K, response_matrix, num_items_per_group = J_g, c = c_init,
#                                    method_lcm = "poLCA", method_dist = "euclidean", #"poLCA"
#                                    method_hclust = "single", method_add_root = "min_cor",
#                                    alpha=0, theta=0, maxiter=100,
#                                    tol=1e-5, na.rm=FALSE, nrep=10, verbose=FALSE, calc.se=TRUE), silent=TRUE)
#   if ('try-error' %in% class(initialization)) initialization_num <- sample.int(9000, size=1)
#   else noError = T
# }

priors <- initialization$priors
initials <- initialization$initials
# plot(initialization$initials$tree_phylo4d, plot.data = F, show.node = T)
# plot(initial_tree, show.node=T, plot.data = F)
# plot(tree, show.node=T, plot.data = F)
cat("c_init = ", initials$c, "\n")
fix_tree_at <- "FALSE"
if (fix_tree_at == "truth"){
  initial_tree <- tree_with_profiles
  # initial_tree@data[as.character(1:K), 1:J] <- logit_response_probs_polca[hungarian_polca,]
  initials$tree_phylo4d <- initial_tree
  fix_tree = TRUE
} else if (fix_tree_at == "misspecified"){
  fix_tree = TRUE
  set.seed(22)
  tree_coal <- as.phylo(rcoal(K, tip.label = paste0("v", 1:K), rooted = TRUE))
  total_length <- nodeHeight(phylo4d(tree_coal), "v1", "root")
  # since the resulting tree is not rooted, we need to sample a divergence time
  # for the root node and add to the tree
  root_edge_length <- 0.2
  # since the total depth rom hclust is not 1, we need to normalize. The tree
  # from hclust will have total depth 1 - root_edge_length
  tree_coal$edge.length <- tree_coal$edge.length / total_length * (1 - root_edge_length)
  # add root node
  tree_coal <- add_root(tree_coal, root_edge_length = root_edge_length, 
                        root_label = "u1", leaf_label = "u2")
  tree_coal <- phylo4d(tree_coal)
  # add node labels
  tipLabels(tree_coal) <- paste0("v", 1:K)
  nodeLabels(tree_coal) <- paste0("u", 1:K)
  # plot(tree_coal, show.node = T)
  tree_coal@data <- initial_tree@data
  initials$tree_phylo4d <- tree_coal
} else if (fix_tree_at == "FALSE"){
  fix_tree = FALSE
}
total_iters <- 8000
burnin <- 5000
controls <- list(
  fix_tree = fix_tree,
  is_linear = is_linear # whether divergence function is linear
)
system.time({
  result <- ddt_lcm_mcmc(K, data = response_matrix, num_items_per_group = J,
                         initials = initials, priors = priors,
                         total_iters = total_iters, controls = controls)
})


### posterior summaries ----------------------------------------------------------

# post hoc label switching
map_index <- which.max(result$loglikelihood_lcm[(burnin+1):total_iters]) + burnin
response_probs_samples <- array(t(result$response_probs_samples[,(burnin+1):total_iters]),
                                dim = c(total_iters-burnin, K, J))
ls_lcm <- label.switching(
  method = c("ECR"),
  zpivot = result$Z_samples[,map_index],
  # mxN integer array of the latent allocation vectors generated from an MCMC algorithm
  z = t(result$Z_samples[,(burnin+1):total_iters]),
  K = K,
  # KxJ array containing the parameter that will be used as a pivot
  prapivot = matrix(result$response_probs_samples[,map_index], nrow=K),
  constraint = 1,
  mcmc = response_probs_samples
)
for (iter in 1:(total_iters-burnin)){
  response_probs_samples[iter,,] <- response_probs_samples[iter, ls_lcm$permutations$`ECR`[iter,],]
}
response_probs_samples_permuted <- response_probs_samples


###  produce summary -----------------------------------
clinear <- TRUE
rmse <- function(x) sqrt(mean(x**2))
l1norm <- function(x) mean(abs(x))
frobenius_norm <- function(x) sqrt(sum(x[lower.tri(x, diag = T)]**2))

s <- summary(result, burnin = burnin, relabel = T) #8000

# initialize data matrix to store the results
dat_result_cols <- c("seed_theta", "seed_Y", "fixtree", "method_hclust",
                     "rmse", "tree_frobenius", "ari", "label", "N")
dat_result <- matrix(0, nrow = 1, ncol = length(dat_result_cols))
colnames(dat_result) <- dat_result_cols
row_idx <- 1
dat_result[row_idx, 1:3] <- c(seed_theta, seed_Y, fix_tree_at)
dat_result[row_idx, c("label", "N")] <- c("DDT-LCM (homogeneous var)", N)

## compare response profile with the truth
# get posterior
response_probs <- matrix(s$response_probs_summary[,"50%"], nrow = K)
# response_probs <-
hungarian_est <- solve_LSAP(sim_data$response_prob %*% t(response_probs), maximum=TRUE)
response_probs_mse <- rmse(response_probs[hungarian_est,] - sim_data$response_prob)
dat_result[row_idx, "rmse"] <- response_probs_mse

# class memberships
class_predictive_probs <- predict_class(s, sim_data$response_matrix)$predictive_probs
class_prediction <- apply(class_predictive_probs, 1, which.max)
dat_result[row_idx, "ari"] <- ARI(class_prediction, sim_data$class_assignment)

# MAP tree
if (fix_tree_at == "FALSE"){
  tree_Sigma_true <- logllk_ddt(1, clinear, Sigma_by_group=rep(1, G), sim_data$tree_with_data, num_items_per_group,
                                tree_structure_old = NULL, dist_mat_old = NULL)$dist_mat$dist_mat_g
  tree_map <- s$tree_map
  maptree_index <- which.max(result$loglikelihood[(burnin+1):total_iters])
  #switch_label
  hungarian_tree <- solve_LSAP(response_probs_samples_permuted[maptree_index, 
                                                               as.numeric(gsub("v", "", tipLabels(tree_map))),] %*% t(sim_data$response_prob), maximum=TRUE)
  # hungarian_tree <- solve_LSAP(sim_data$response_prob %*% t(response_probs_samples_permuted[maptree_index, as.numeric(gsub("v", "", tipLabels(tree_map))),] ), maximum=TRUE)
  tipLabels(tree_map) <- paste0("v", hungarian_tree)
  # plot(tree_map, show.node = T, plot.data = F)
  # plot(tree_with_profiles, show.node = T, plot.data = F)
  tree_Sigma_est <- logllk_ddt(1, clinear, Sigma_by_group, tree_map, num_items_per_group,
                               tree_structure_old = NULL, dist_mat_old = NULL)$dist_mat$dist_mat_g
  dat_result[row_idx, "tree_frobenius"] <- frobenius_norm(tree_Sigma_est-tree_Sigma_true)
} else{
  dat_result[row_idx, "tree_frobenius"] <- NA
}
summarized_result_ddtlcm_homovar <- dat_result

save_dir <- "simulation_semisynthetic/data/1_semi_result"
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/N", N, "_J", J, "_K", K)
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/seedtheta", seed_theta)
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/ddtlcm_homovar")
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
file_name <- paste0(save_dir, "/simdat_seed", seed_Y, ".RData")
save(seed_theta, seed_Y, sim_data, initialization, initials, 
     result, summarized_result_ddtlcm_homovar, file = file_name)

