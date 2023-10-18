test_that("'H_n' returns valid Harmonic series.", {
  expect_equal(H_n(0), 0)
  expect_equal(H_n(1), 1)
  expect_equal(H_n(2), 1.5)
  expect_equal(H_n(4), 1+1/2+1/3+1/4)
  expect_equal(H_n(-1), Inf)
  expect_error(H_n(Inf))
  expect_error(H_n(c(1,4,5)))
})


