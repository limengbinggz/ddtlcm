test_that("'summary.ddt_lcm' returns a valid summary of DDT-LCM", {
  # load the result of fitting semi-synthetic data with 100 (for the sake of time) posterior samples
  data(result_diet_1000iters)
  summarized_result <- summary(result_diet_1000iters, burnin = 50, relabel = TRUE, be_quiet = TRUE)
  K <- result_diet_1000iters$setting$K
  J <- length(unlist(result_diet_1000iters$setting$item_membership_list))
  N <- nrow(result_diet_1000iters$data)

  expect_length(summarized_result, 10)
  expect_s3_class(summarized_result, "summary.ddt_lcm")
  expect_s4_class(summarized_result$tree_map, "phylo4d")
  expect_equal(summarized_result$setting$K, K)
  expect_equal(nrow(summarized_result$response_probs_summary), K*J)

})
