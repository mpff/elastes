data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))
data_curve2 <- data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
data_curves <- list(data_curve1, data_curve2)


test_that("Test computed variance and distances in [0,1]", {
  mean <- compute_elastic_shape_mean(data_curves)
  expect_gt(mean$variance, 0)
  expect_lte(mean$variance, 1)
  for (d in mean$distances) {
    expect_gte(d, 0)
    expect_lte(d, 1)
  }
})


test_that("Test variance == sum-of-squared distances / N", {
  mean <- compute_elastic_shape_mean(data_curves, max_iter = 0)
  ssdn <- sum(mean$distances^2)*2/(length(mean$distances)*(length(mean$distances) - 1))
  expect_equal(mean$variance, ssdn, tolerance = 1e-1)
})


test_that("Test unelastic full procrustes mean has distances", {
  mean <- compute_elastic_shape_mean(data_curves, max_iter = 0)
  expect_equal(mean$data_curves[[1]]$t, mean$data_curves[[1]]$t_optim)
  expect_false(is.null(attr(mean$data_curves[[1]], "dist_to_mean")))
  mean2 <- compute_elastic_shape_mean(data_curves, max_iter = 0, pfit_method = "polygon")
  expect_equal(mean2$data_curves[[1]]$t, mean2$data_curves[[1]]$t_optim)
  expect_false(is.null(attr(mean2$data_curves[[1]], "dist_to_mean")))
})


test_that("Initial rot, scaling do not matter much", {
  knots <- seq(0, 1, length = 11)

  # Apply rot, scaling to data_curve2
  beta <- 2.5
  theta <- 0.7*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  data_curves <- list(data_curve1, data_curve2)

  mean1 <- compute_elastic_shape_mean(data_curves, knots = knots, type = "smooth", penalty = 2)
  mean2 <- compute_elastic_shape_mean(data_curves, knots = knots, type = "polygon", penalty = 2)

  # Apply different rot, scaling to data_curve2
  beta <- 0.25
  theta <- 2.3*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  data_curves <- list(data_curve1, data_curve2)

  # Get means "smooth"
  mean3 <- compute_elastic_shape_mean(data_curves, knots = knots, type = "smooth", penalty = 2)
  mean4 <- compute_elastic_shape_mean(data_curves, knots = knots, type = "polygon", penalty = 2)

  expect_equal(mean1$coefs, mean3$coefs, tolerance=1e-1)
  expect_equal(mean2$coefs, mean4$coefs, tolerance=1e-1)
})
