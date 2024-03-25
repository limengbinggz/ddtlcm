test_that("'plot' returns a valid plot object.", {
  data(result_diet_1000iters)
  burnin <- 500

  # trace plots
  expect_error(plot(x = result_diet_1000iters, c("resonseprob_1,1,1"), burnin = 500),
               "parameter_names[1] invalid. Please check the appropriate format in the function description", fixed=TRUE)
  expect_error(plot(x = result_diet_1000iters, c("classsprob_1,1,1"), burnin = 500),
               "parameter_names[1] invalid. Please check the appropriate format in the function description", fixed=TRUE)
  expect_error(plot(x = result_diet_1000iters, c("cc"), burnin = 500),
               "parameter_names[1] invalid. Please check the appropriate format in the function description", fixed=TRUE)
  expect_error(plot(x = result_diet_1000iters, c("1_diffusionvar"), burnin = 500),
               "parameter_names[1] invalid. Please check the appropriate format in the function description", fixed=TRUE)
  expect_error(plot(x = result_diet_1000iters, c("responseprob_1,1,1"), burnin = 1001),
               "The number of burn-in posterior samples, burnin, should be a positive integer smaller than the total number of samples.", fixed=TRUE)

  summarized_result <- summary(result_diet_1000iters, burnin, relabel = TRUE, be_quiet = TRUE)
  plot_out1 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "all")
  expect_s3_class(plot_out1, c("gg", "ggplot", "ggarrange"))

  plot_out2 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "tree")
  expect_s3_class(plot_out2, c("gg", "ggplot", "ggtree"))

  plot_out3 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "profile")
  expect_s3_class(plot_out3, c("gg", "ggplot"))
})


