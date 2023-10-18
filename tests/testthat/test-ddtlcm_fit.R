test_that("'ddtlcm_fit' returns a valid 'ddt_lcm' object.", {
  
  # load the MAP tree structure obtained from the real HCHS/SOL data
  data(data_synthetic)
  # extract elements into the global environment
  list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
  # run DDT-LCM
  total_iters = 5
  result <- ddtlcm_fit(K = 3, data = response_matrix, item_membership_list, total_iters = total_iters)
  K <- result$setting$K
  J <- length(unlist(result$setting$item_membership_list))
  N <- nrow(result$data)
  G <- length(item_membership_list)
  
  expect_s3_class(result, "ddt_lcm")
  expect_equal(length(result), 11L)
  expect_equal(dim(result$response_probs_samples), c(total_iters, K, J))
  expect_equal(dim(result$class_probs_samples), c(K, total_iters))
  expect_equal(dim(result$Z_samples), c(N, total_iters))
  expect_equal(dim(result$Sigma_by_group_samples), c(G, total_iters))
})


