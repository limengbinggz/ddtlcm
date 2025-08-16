#######################################################
#' Simulate data along different DDT trees with multiple
#' divergence parameters, from early to late divergence
#' Generate different sets of theta along the given tree.
#' For each set of theta, generate multiple multivariate
#' binary responses.
#' Independent normal priors on the logistic-transformed
#' item responsibilities. Group-specific variance parameters.
#' Z_i ~ Cat(\pi)
#' P(Y_ij = 1 \mid Z_i = k) = \theta_{kj}
#' logistic(\theta_{kj}) ~ N(0, \sigma^2), for all item j.
#' #######################################################
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
#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

tree_idx = as.integer(args[1])
N <- as.integer(args[2])
seed_theta = as.integer(args[3])
seed_Y = as.integer(args[4])

## specify parameters -----------------------------
tree_idx <- 3
seed_theta <- 1
seed_Y <- 2
N=400

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
cat("c_init = ", initials$c, "\n")
leaf_data <- as.matrix(initials$tree_phylo4d@data[as.character(1:K), 1:J])
initials$item_probability <- expit(leaf_data)
total_iters <- 80
burnin <- 50

source("functions/Gibbs_bayeslcm.R")
controls <- list()
system.time({
  result <- gibbs_bayeslcm(K, data = response_matrix, num_items_per_group = J,
                         initials = initials, priors = priors,
                         total_iters = total_iters, controls = controls)
})


###  produce summary -----------------------------------
clinear <- TRUE
rmse <- function(x) sqrt(mean(x**2))
l1norm <- function(x) mean(abs(x))
frobenius_norm <- function(x) sqrt(sum(x[lower.tri(x, diag = T)]**2))
tree_Sigma_true <- logllk_ddt(1, clinear, Sigma_by_group=Sigma_by_group, sim_data$tree_with_data, num_items_per_group,
                              tree_structure_old = NULL, dist_mat_old = NULL)$dist_mat$dist_mat_g

# linkage methods to perform post hoc agglomerative clustering
method_hclust_list <- c("ward.D", "ward.D2", "single", "average", "complete", "mcquitty")
n_method_hclust <- length(method_hclust_list)

## summarize result
tryCatch({
  s <- summary(result, burnin = burnin, relabel = T) 
}, error=function(cond) {
  s <- summary(result, burnin = burnin, relabel = F) 
})

fix_tree_at <- "FALSE"
# initialize data matrix to store the results
dat_result_cols <- c("seed_theta", "seed_Y", "fixtree", "method_hclust",
                     "rmse", "tree_frobenius", "ari", "label", "N")
dat_result <- matrix(0, nrow = n_method_hclust, ncol = length(dat_result_cols))
colnames(dat_result) <- dat_result_cols
row_idx <- 1
dat_result[row_idx, 1:3] <- c(seed_theta, seed_Y, fix_tree_at)
dat_result[row_idx, c("label", "N")] <- c("BayesLCM (homogeneous var) + HC", N)

## compare response profile with the truth
# get posterior
response_probs <- matrix(s$response_probs_summary[,"50%"], nrow = K)
hungarian_est <- solve_LSAP(sim_data$response_prob %*% t(response_probs), maximum=TRUE)
response_probs_mse <- rmse(response_probs[hungarian_est,] - sim_data$response_prob)
dat_result[row_idx, "rmse"] <- response_probs_mse

# class memberships
class_predictive_probs <- predict_class(s, sim_data$response_matrix)$predictive_probs
class_prediction <- apply(class_predictive_probs, 1, which.max)
ari <- ARI(class_prediction, sim_data$class_assignment)
dat_result[row_idx, "ari"] <- ari

## compare estimated tree with the truth
# need to re-create the tree here to remove the initial branch
# create distance matrix between classes
response_probs <- response_probs[hungarian_est,]
rownames(response_probs) <- paste0("v", 1:K)
# create distance matrix between classes
leaf_dist <- dist(logit(response_probs), method = "euclidean")
# hierarchical clustering
for (method_hclust in method_hclust_list){
  hclust_phylo <- as.phylo(hclust(leaf_dist, method = method_hclust))
  total_length <- nodeHeight(phylo4d(hclust_phylo), "v1", "root")
  hclust_phylo$edge.length <- hclust_phylo$edge.length / total_length
  tree_phylo4d <- phylo4d(hclust_phylo)
  tipLabels(tree_phylo4d) <- paste0("v", 1:K)
  # switch labels for trees
  hungarian_tree <- solve_LSAP(
    response_probs[as.numeric(gsub("v", "", tipLabels(tree_phylo4d))),] %*% t(sim_data$response_prob), maximum=TRUE)
  tipLabels(tree_phylo4d) <- paste0("v", hungarian_tree)
  
  nodeLabels(tree_phylo4d) <- paste0("u", 1:(K-1))
  tree_Sigma_hclust <- create_leaf_cor_matrix(tree_phylo4d, mc.cores = 1L)
  if (method_hclust == method_hclust_list[1]){
    dat_result[row_idx, "method_hclust"] <- method_hclust
    dat_result[row_idx, "tree_frobenius"] <- frobenius_norm(tree_Sigma_hclust - tree_Sigma_true)
  } else{
    frob <- frobenius_norm(tree_Sigma_hclust - tree_Sigma_true)
    dat_result[row_idx, ] <- c(seed_theta, seed_Y, fix_tree_at, 
                               method_hclust, response_probs_mse, frob, ari,
                               "BayesLCM (homogeneous var) + HC", N)
  }
  row_idx <- row_idx + 1
}

summarized_result_bayeslcm_homovar <- dat_result

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
save_dir <- paste0(save_dir, "/bayeslcm_homovar")
if (!dir.exists(save_dir)){teee9
  print("create")
  dir.create(save_dir)
}
file_name <- paste0(save_dir, "/simdat_seed", seed_Y, ".RData")
save(sim_data, seed_theta, seed_Y, result, 
     tree_idx, summarized_result_bayeslcm_homoovar, Sigma_dat, file = file_name)


