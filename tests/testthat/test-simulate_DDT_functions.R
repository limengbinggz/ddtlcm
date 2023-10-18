test_that("'simulate_DDT_tree' returns a valid simulated 'phylo' tree object", {
  options(warn = -1)
  set.seed(6)
  K <- 6
  c <- 5
  c_order <- 1
  tree1 <- simulate_DDT_tree(K, c, c_order)
  expect_s3_class(tree1, 'phylo')
  expect_equal(tree1$Nnode, K)
  expect_equal(tree1$tip.label, paste0("v", 1:K))

  # check leaf divergence time
  tree1_phylo4 <- as(tree1, "phylo4")
  expect_equal(nodeHeight(tree1_phylo4, "v1", "root"), 1)
  
  expect_s3_class(simulate_DDT_tree(K, c, c_order, alpha = 0.4, theta = 0.1), 'phylo')
})


test_that("'simulate_lcm_response' returns a valid matrix of binary response matrix", {
  set.seed(6)
  # number of latent classes
  K <- 6
  # number of items
  J <- 78
  response_prob <- matrix(runif(K*J), nrow = K)
  class_probability <- rep(1/K, K)
  # number of individuals
  N <- 100
  response_matrix <- simulate_lcm_response(N, response_prob, class_probability)
  
  expect_length(response_matrix, 2)
  expect_equal(unique(c(response_matrix$response_matrix)), 0:1)
  expect_equal(sort(unique(c(response_matrix$class_assignment))), 1:K)
})


test_that("'simulate_lcm_given_tree' returns a valid list of simulated binary response matrix and parameters", {
  set.seed(6)
  # load the MAP tree structure obtained from the real HCHS/SOL data
  data(parameter_diet)
  # unlist the elements into variables in the global environment
  list2env(setNames(parameter_diet, names(parameter_diet)), envir = globalenv())
  # number of individuals
  N <- 496
  # set random seed to generate node parameters given the tree
  seed_parameter = 1
  # set random seed to generate multivariate binary observations from LCM
  seed_response = 1
  # simulate data given the parameters
  sim_data <- simulate_lcm_given_tree(tree_phylo, N,
                                      class_probability, item_membership_list, Sigma_by_group,
                                      root_node_location = 0, seed_parameter = 1, seed_response = 1)
  expect_length(sim_data, 7)
  expect_s4_class(sim_data$tree_with_parameter, 'phylo4d')
})

