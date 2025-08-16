#######################################################
#' Simulate data along the MAP tree (K=6) estimated from real
#' data. 
#######################################################

curr_dir <- "/Users/mengbing/Dropbox (University of Michigan)/from_box/research/tree_multivariate_binary/reproducibility_codes"
setwd(curr_dir)

library(ggtree)
source("functions/utils.R")
source("functions/simulate_DDT_functions.R")
options(warn = -1)

load("hchs_simulated_data/tree_real_nodata.RData") 
# sample size in the HCHS/SOL study South American subpopulation
N = 496
# number of classes
K <- 6

### extract parameters -----------------------------------------
num_items_per_group <- J_g
G <- length(J_g)
item_group_membership <- rep(1:G, J_g)
# total number of items
J <- sum(J_g) #78
class_probability_true <- class_probability

### Simulate data -----------------------------------------------------------

## first simulate class profiles along the given tree, following Gaussian diffusion processes
# location of root node in the DDT
root_node_location0 <- rep(0, J)
# set seed for generating class profiles
seed_theta = 1
# is_linear <- ifelse(is_linear == "TRUE", T, F)
set.seed(seed_theta)
tree_phylo <- as.phylo(extractTree(tree))
as.phylo(extractTree(tree))
tree_with_profiles <- gen_location_with_selector(tree_phylo, Sigma_by_group, item_group_membership, 
                                                 root_node_location = root_node_location0)

## second simulate multivariate binary responses following latent class model
# extract logit of response profiles by class
leaf_data <- as.matrix(tree_with_profiles@data[as.character(1:K), 1:J])
response_prob <- expit(leaf_data)
# set seed for generating responses
seed_Y = 1
set.seed(seed_Y)
# simulate multivariate responses
sim_responses <- gen_response(N, response_prob, class_probability)
# N x J binary response matrix
response_matrix <- sim_responses$response_matrix
food_groups <- c("diary", "fat", "fruit", "grain", "meat", "sugar", "veg")
item_names <- c()
for (g in 1:G) {
  item_names <- c(item_names, paste0(food_groups[g], "_", 1:J_g[g]))
}
colnames(response_matrix) <- item_names

# save simulated data as a list
sim_data <- list()
sim_data[["tree_with_data"]] <- tree_with_profiles
sim_data[["class_probability"]] <- class_probability
sim_data[["response_prob"]] <- response_prob
sim_data[["response_matrix"]] <- response_matrix
sim_data[["class_assignment"]] <-  sim_responses$class_assignment
sim_data[["Sigma_by_group"]] <- Sigma_by_group

file_name <- paste0("hchs_simulated_data/hchs_simulated_data.RData")
save(seed_theta, seed_Y, num_items_per_group, N, sim_data, file = file_name)



