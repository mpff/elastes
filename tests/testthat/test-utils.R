data_curve <- data.frame(x1 = c(0,0.5,1), x2 = c(0,0,0))

test_that("Test can calculate polygon length", {
  expect_equal(get_polygon_length(data_curve), 1)
})

test_that("Test can calculate curve center", {
  expect_equal(get_center(data_curve),  c(x1 = 0.5, x2 = 0))
})

test_that("Test distance between same curves is zero",{
  srv_data_curve <- elasdics::get_srv_from_points(data_curve)
  knots <- c(0,0.5,1)
  type <- "polygon"
  coefs <- as.matrix(srv_data_curve[,-1])
  distance <- get_distance(
    srv_curve = function(t) t(make_design(t, knots=knots, type=type) %*% coefs),
    s = c(srv_data_curve$t, 1),
    q = t(srv_data_curve[,-1])
  )
  expect_equal(distance, 0, tolerance=1e-6)
})
