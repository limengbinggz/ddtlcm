test_that("'predict' returns a valid prediction.", {
  data(result_diet_shortchain)
  N <- nrow(result_diet_shortchain$data)
  K <- result_diet_shortchain$setting$K
  burnin <- 50
  summarized_result <- summary(result_diet_shortchain, burnin, relabel = TRUE, be_quiet = TRUE)
  predicted <- predict(summarized_result, result_diet_shortchain$data)
  expect_length(predicted, 2)
  expect_length(predicted$class_assignments, N)
  expect_equal(dim(predicted$predictive_probs), c(N, K))

  expect_error(predict(result_diet_shortchain, result_diet_shortchain$data))
})


