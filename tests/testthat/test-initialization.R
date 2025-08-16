test_that("'initialize' returns a valid initialization for DDT-LCM.", {
  set.seed(60)
  data(data_synthetic)
  # extract elements into the global environment
  list2env(setNames(data_synthetic, names(data_synthetic)), envir = globalenv())
  K <- 3
  G <- length(item_membership_list)
  initials <- initialize(K, data = response_matrix, item_membership_list,
                         c=1, c_order=1,  method_lcm = "random",
                         fixed_initials = NULL, fixed_priors = NULL)
  expect_length(initials, 3)
  expect_length(initials$initials, 6)
  expect_length(initials$priors, 6)
  expect_length(initials$model_lcm, 3)
  
  
  fixed_initials <- list("c" = 5)
  fixed_priors <- list("rate_sigma" = rep(3, G), "shape_c" = 2, "rate_c" = 2)
  initials <- initialize(K, data = response_matrix, item_membership_list,
    c=1, c_order=1,  method_lcm = "random",
    fixed_initials = fixed_initials, fixed_priors = fixed_priors)
  expect_equal(initials$initials$c, 5)
  expect_equal(initials$priors$rate_sigma, rep(3, G))
  expect_equal(initials$priors$shape_c, 2)
  expect_equal(initials$priors$rate_c, 2)

  expect_error(initialize(K, data = response_matrix, item_membership_list,
                         c=1, c_order=1,  method_lcm = "other",
                         fixed_initials = fixed_initials, fixed_priors = fixed_priors))
  expect_error(initialize(K, data = response_matrix, item_membership_list,
                         c=1, c_order=1,   method_add_root = NA, 
                         fixed_initials = fixed_initials, fixed_priors = fixed_priors))
})


