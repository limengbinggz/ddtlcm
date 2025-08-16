###########################################################
#' Set up 4 different tree structures with 3 leaves, representing different 
#' degrees of class separation. For each tree, simulate theta for K=3
#' class profiles. For each set of theta, generate multiple multivariate
#' binary responses.
###########################################################
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

args = commandArgs(trailingOnly=TRUE)
tree_idx = as.integer(args[1]) # index of tree structure
N <- as.integer(args[2]) # sample size
seed_theta = as.integer(args[3]) # random seed for theta
seed_Y = as.integer(args[4]) # random seed for data Y


# ## specify parameters -----------------------------
# tree_idx <- 3
# N <- 400
# seed_theta <- 1
# seed_Y <- 1
# fix_tree_at <- "FALSE"


# the number of classes K
K = 3
class_probability <- c(0.4, 0.3, 0.3)
num_items_per_group <- c(rep(10, 5), 15, 15)
# number of groups
G <- length(num_items_per_group)
item_group_membership <- rep(1:G, num_items_per_group)
# total number of items
J <- sum(num_items_per_group) #80
# variance of logit response probabilities of items in each group
Sigma_by_group <- c(rep(0.6**2, 5), rep(2**2, 2))

### Simulate data -----------------------------------------------------------
# location of root node in the DDT
root_node_location0 <- rep(0, J)
is_linear <- T
if (tree_idx == 1){
  tr_txt <- "(((v1:0.85, v2:0.85):0.1, v3:0.95):0.05);"
  tree <- read.tree(text = tr_txt)
  tree$node.label <- paste0("u", 1:Nnode(tree))
} else if (tree_idx == 2){
  tr_txt <- "(((v1:0.5, v2:0.5):0.25, v3:0.75):0.25);"
  tree <- read.tree(text = tr_txt)
  tree$node.label <- paste0("u", 1:Nnode(tree))
  # plot(tree, show.node.label = T)
} else if (tree_idx == 3){
  tr_txt <- "(((v1:0.25, v2:0.25):0.65, v3:0.9):0.1);"
  tree <- read.tree(text = tr_txt)
  tree$node.label <- paste0("u", 1:Nnode(tree))
  plot(tree, show.node.label = T)
} else if (tree_idx == 4){
  tr_txt <- "(((v1:0.25, v2:0.25):0.15, v3:0.4):0.6);"
  tree <- read.tree(text = tr_txt)
  tree$node.label <- paste0("u", 1:Nnode(tree))
  # plot(tree, show.node.label = T)
}

set.seed(seed_theta)
tree_with_profiles <- gen_location_with_selector(tree, Sigma_by_group, item_group_membership,
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
initialization <- initialize(K, response_matrix, num_items_per_group = J, c = c_init,
                             method_lcm = "random", method_dist = "euclidean", 
                             method_hclust = "single", method_add_root = "min_cor")

# # alternatively, we may initialize using poLCA
# # poLCA may fail. If this is the case, use another seed
# noError = F
# initialization_num <- 1
# repeat{
#   if (noError == T) break
#   set.seed(initialization_num)
#   initialization <- try(initialize(K, response_matrix, num_items_per_group = num_items_per_group, c = c_init,
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
cat("c_init = ", initials$c, "\n")
fix_tree = FALSE
controls <- list(
  fix_tree = fix_tree,
  is_linear = is_linear
)

total_iters <- 80
burnin <- 50
system.time({
  result <- ddt_lcm_mcmc(K, data = response_matrix, num_items_per_group = J,
                         initials = initials, priors = priors,
                         total_iters = total_iters, controls = controls)
})



### posterior summaries ----------------------------------------------------------
map_index <- which.max(result$loglikelihood_lcm[(burnin+1):total_iters]) + burnin
# switch labels. m = # of iterations
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

s <- summary(result, burnin = burnin, relabel = T) 

# initialize data matrix to store the results
dat_result_cols <- c("seed_theta", "seed_Y", "tree", "fixtree", "method_hclust",
                     "rmse", "tree_frobenius", "ari", "label", "N")
dat_result <- matrix(0, nrow = 1, ncol = length(dat_result_cols))
colnames(dat_result) <- dat_result_cols
row_idx <- 1
dat_result[row_idx, 1:4] <- c(seed_theta, seed_Y, tree_idx, fix_tree_at)
dat_result[row_idx, c("label", "N")] <- c("DDT-LCM (homogeneous var)", N)

## compare response profile with the truth
# get posterior
response_probs <- matrix(s$response_probs_summary[,"50%"], nrow = K)
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
  tipLabels(tree_map) <- paste0("v", hungarian_tree)
  tree_Sigma_est <- logllk_ddt(1, clinear, Sigma_by_group, tree_map, num_items_per_group,
                               tree_structure_old = NULL, dist_mat_old = NULL)$dist_mat$dist_mat_g
  dat_result[row_idx, "tree_frobenius"] <- frobenius_norm(tree_Sigma_est-tree_Sigma_true)
} else{
  dat_result[row_idx, "tree_frobenius"] <- NA
}

summarized_result_ddtlcm_homovar <- dat_result


### Summarize variance parameter estimation ---------------------
Sigma_by_group_true <- sim_data$Sigma_by_group
bias <- s$Sigma_summary[, "Mean"] - Sigma_by_group_true
cover <- Sigma_by_group_true >= s$Sigma_summary[, "2.5%"] & Sigma_by_group_true <= s$Sigma_summary[, "97.5%"] 
Sigma_dat <- data.table(tree_idx = tree_idx, N = N, seed_theta = seed_theta, seed_Y = seed_Y,
                        var_name = names(bias), 
                        est = s$Sigma_summary[, "Mean"], sd = s$Sigma_summary[, "SD"], 
                        bias = bias, cover = cover)



save_dir <- "simulation_synthetic/data/1_syn_result"
if (!dir.exists(save_dir)){
  print("create")
  dir.create(save_dir)
}
save_dir <- paste0(save_dir, "/tree", tree_idx, "_N", N, "_J", J, "_K", K)
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
save(sim_data, seed_theta, seed_Y, result, 
     tree_idx, summarized_result_ddtlcm_homovar, Sigma_dat, file = file_name)
