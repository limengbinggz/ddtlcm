test_that("'predict' returns a valid prediction.", {
  data(result_diet)
  N <- nrow(result_diet$data)
  K <- result_diet$setting$K
  burnin <- 50
  summarized_result <- summary(result_diet, burnin, relabel = TRUE, be_quiet = TRUE)
  predicted <- predict(summarized_result, result_diet$data)
  expect_length(predicted, 2)
  expect_length(predicted$class_assignments, N)
  expect_equal(dim(predicted$predictive_probs), c(N, K))
  
  expect_error(predict(result_diet, result_diet$data))
})


