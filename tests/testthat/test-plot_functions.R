test_that("'plot' returns a valid plot object.", {
  data(result_diet_1000iters)
  burnin <- 50
  summarized_result <- summary(result_diet_1000iters, burnin, relabel = TRUE, be_quiet = TRUE)
  plot_out1 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "all")
  expect_s3_class(plot_out1, c("gg", "ggplot", "ggarrange"))

  plot_out2 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "tree")
  expect_s3_class(plot_out2, c("gg", "ggplot", "ggtree"))

  plot_out3 <- plot(x = summarized_result, item_name_list = NULL, plot_option = "profile")
  expect_s3_class(plot_out3, c("gg", "ggplot"))
})


