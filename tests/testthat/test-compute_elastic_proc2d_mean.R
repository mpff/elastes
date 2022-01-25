
test_that("Test type = 'smooth': Initial rot, scaling do not matter much", {
  knots <- seq(0, 1, length = 11)
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))  # unscaled, unrotated

  # Apply rot, scaling to data_curve2
  beta <- 2.5
  theta <- 0.7*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean1 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "smooth", penalty = 2)
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }

  # Apply different rot, scaling to data_curve2
  beta <- 0.25
  theta <- 2.3*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean3 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "smooth", penalty = 2)
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }
  expect_equal(mean1$coefs, mean3$coefs, tolerance=1e-1)
})


test_that("Test type = 'polygon': initial rot, scaling do not matter much", {
  knots <- seq(0, 1, length = 16)
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))  # unscaled, unrotated

  # Apply rot, scaling to data_curve2
  beta <- 2.5
  theta <- 0.7*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean2 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "polygon", penalty = 2)
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }
  # Apply different rot, scaling to data_curve2
  beta <- 0.25
  theta <- 2.3*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean4 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "polygon", penalty = 2)
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }
  expect_equal(mean2$coefs, mean4$coefs, tolerance=1e-1)
})


test_that("Test pfit_method = 'polygon' (depriciated): initial rot, scaling do not matter much", {
  knots <- seq(0, 1, length = 11)
  data_curve1 <- data.frame(x1 = sin(1:7/4*pi), x2 = cos(1:7/4*pi))  # unscaled, unrotated

  # Apply rot, scaling to data_curve2
  beta <- 2.5
  theta <- 0.7*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean1 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "smooth", penalty = 2, pfit_method="polygon")
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }

  # Apply different rot, scaling to data_curve2
  beta <- 0.25
  theta <- 2.3*pi
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2, ncol = 2)
  data_curve2 <- data_curve <- beta * data.frame(x1 = sin(1:15/8*pi), x2 = cos(1:15/8*pi))
  data_curve2 <- as.matrix(data_curve2) %*% t(rot)
  data_curve2 <- as.data.frame(data_curve2)
  colnames(data_curve2) <- c("x1","x2")
  # Get means
  data_curves <- list(data_curve1, data_curve2)
  w <- capture_warnings(
    mean3 <- compute_elastic_proc2d_mean(data_curves, knots = knots, type = "smooth", penalty = 2, pfit_method="polygon")
  )
  if (length(w) > 0) {
    for (warn in w) expect_equal(warn, "there is *no* information about some basis coefficients")
  }
  expect_equal(mean1$coefs, mean3$coefs, tolerance=1e-1)
})
