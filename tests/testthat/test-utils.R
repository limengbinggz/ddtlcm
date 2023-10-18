test_that("'a_t_one' corrrectly computes divergence function a(t) = c / (1-t).", {
  expect_equal(a_t_one(1, 0.5), 2)
  expect_equal(a_t_one(5, 0.8), 25)
})

test_that("'a_t_one_cum' corrrectly computes cumulative hazard function for a(t) = c / (1-t).", {
  expect_equal(a_t_one_cum(1, 0.99), -log(0.01))
  expect_equal(a_t_one_cum(2, 0.99), -2*log(0.01))
})

test_that("'A_t_inv_one' corrrectly computes cumulative hazard function for a(t) = c / (1-t).", {
  expect_equal(A_t_inv_one(1, 2), 1 - exp(- 2))
})

test_that("expit function corrrect.", {
  expect_equal(expit(0), 0.5)
  expect_equal(expit(-Inf), 0)
  expect_equal(expit(Inf), 1)
  expect_equal(expit(c(-1, 1)), 1 / (1 + exp(-c(-1,1))))
})

test_that("logit function corrrect.", {
  expect_equal(logit(0.5), 0)
  expect_equal(logit(1), 5)
  expect_equal(logit(0), -5)
  expect_equal(logit(c(0.1)),  log(0.1 / 0.9))
})

test_that("'create_leaf_cor_matrix' returns a valid tree-structured covariance matrix.", {
  tr_txt <- "(((v1:0.85, v2:0.85):0.1, v3:0.95):0.05);"
  tree <- read.tree(text = tr_txt)
  tree$node.label <- paste0("u", 1:Nnode(tree))
  tree_phylo4d <- as(tree, "phylo4d")
  correct_matrix <- matrix(c(1, 0.15, 0.05, 
                             0.15, 1, 0.05, 
                             0.05, 0.05, 1.0), nrow = 3, byrow = TRUE)
  rownames(correct_matrix) <- paste0("v", 1:3)
  colnames(correct_matrix) <- paste0("v", 1:3)
  expect_equal(create_leaf_cor_matrix(tree_phylo4d), correct_matrix)
})












