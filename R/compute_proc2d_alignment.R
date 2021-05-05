#' @title Optimal rotation and scaling alignment to a smooth curve
#' @description Finds optimal rotation and scaling alignment for a discrete open srv curve to a smooth curve
#' @param q complex srv curve with parametrization, needs to be vectorised.
#' The result of a call to \code{get_model_data_complex}
#' @param coefs.compl complex coefficients of smooth curve
#' @param type spline degree
#' @param knots basis knots
#' @return optimal rotation G and scaling b
#' @importFrom stats approx

compute_proc2d_alignment <- function(q, coefs.compl, type, knots, h){

  arg.grid = seq(0, 1, by=h)

  # Calculate overlap of arg.grid and t_optims.
  idx <- findInterval(arg.grid, q$m_long)
  idx.bool <- which(idx > 0 & idx < length(q$m_long))
  arg.grid.x <- arg.grid[idx.bool]

  # Linear interpolation of srv data curve on overlap.
  q_approx_x <- approx(x=q$m_long, y=Re(q$q_m_long), xout=arg.grid.x)$y
  q_approx_y <- approx(x=q$m_long, y=Im(q$q_m_long), xout=arg.grid.x)$y
  q_approx <- complex(real=q_approx_x, imaginary=q_approx_y)

  # Evaluate mean function on overlap.
  mean_eval <- make_design(arg.grid.x, knots = knots, type = type) %*% coefs.compl

  # Calculate pfit scaling+rotation
  qm <- Conj(q_approx) * mean_eval
  qq <- Conj(q_approx) * q_approx

  # Numerical Integration using trapezoid rule.
  qm <- h * ( sum(qm) - 0.5*qm[1] - 0.5*qm[length(qm)] )
  qq <- h * ( sum(qq) - 0.5*qq[1] - 0.5*qq[length(qq)] )

  # Calculate G and b
  list("qm" = qm, "qq" = qq, "pfit_coefs" = NULL)
}
