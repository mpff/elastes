knots <- c(0, 1)

test_that("Test can get inner and outer knots for both types", {
  expect_equal(get_knots(knots, "polygon"), c(0, 1))
  expect_equal(get_knots(knots, "smooth"), c(-1, 0, 1, 2))
})

test_that("Test can get design matrix for both types", {
  t <- c(0, 1)
  dp <- matrix(c(1, 1), nrow=2)
  expect_equal(make_design(t, knots, "polygon"), dp)
  ds <- matrix(c(1, 0, 0, 1), nrow=2)
  expect_equal(make_design(t, knots, "smooth"), ds)
})

test_that("Test can get Gram matrix for equidistant knots", {
  gp <- matrix(1)
  expect_equal(get_gram_matrix(knots, "polygon"), gp)
  gs <- matrix(c(1/3, 1/6, 1/6, 1/3), nrow=2)
  expect_equal(get_gram_matrix(knots, "smooth"), gs)
})

test_that("Test can get Gram matrix for type 'polygon' with irregular knots", {
  knots_irreg <- c(0, 0.1, 1)
  gip <- matrix(c(0.1, 0, 0, 0.9), nrow=2)
  expect_equal(get_gram_matrix(knots_irreg, "polygon"), gip)
})
