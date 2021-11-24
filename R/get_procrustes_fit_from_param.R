#' Helper functions for calculating procrustes data curve from rotation, scaling and translation parameters.
#' @name get_procrustes_fit_from_param
#' @description Compute the procrustes fit given optimal rotation, scaling and translation.
#' @param data_curve A \code{data.frame} with observed points on a curve.
#' Each row is one point, each variable one coordinate direction. If there is a variable \code{t},
#' it is treated as the time parametrisation, not as an additional coordinate.
#' @param rot The rotation (in radian).
#' @param scale The scaling.
#' @param trans The translation.
#' @export

get_procrustes_fit_from_param <- function(data_curve, rot, scale, trans)
{
  names <- colnames(data_curve)
  if("t" %in% names) {
    t <- data_curve$t
    data_curve <- data_curve[, names(data_curve) != "t"]
  }
  if("t_optim" %in% names) {
    t_optim <- data_curve$t
    data_curve <- data_curve[, names(data_curve) != "t_optim"]
  }

  # Remove translation.
  data_curve <- data_curve - matrix(trans, nrow = nrow(data_curve),
                                    ncol = ncol(data_curve), byrow = TRUE)
  # Remove rotation.
  mat <- matrix(c(cos(-rot), - sin(-rot), sin(-rot), cos(-rot)), nrow = 2, ncol = 2)
  data_curve.rot <- as.matrix(data_curve) %*% t(mat)
  # Remove scaling.
  data_curve.rot.scale <- 1/scale * data_curve.rot

  data_curve <- as.data.frame(data_curve.rot.scale)
  colnames(data_curve) <- names[!names %in% c("t", "t_optim")]
  if("t" %in% names) {
    data_curve$t <- t
  }
  if("t_optim" %in% names) {
    data_curve$t_optim <- t_optim
  }
  data_curve <- data_curve[,names]
  return(data_curve)
}

#' Helper functions for calculating procrustes data curve from rotation, scaling and translation parameters.
#' @name get_procrustes_fit
#' @description Compute the procrustes fit.
#' @param data_curve A \code{data.frame} in an \code{elastic_proc2d_mean} object.
#' @export

get_procrustes_fit <- function(data_curve)
{
  # Get paramaters.
  trans <- attr(data_curve, "translation")
  rot <- attr(data_curve, "rotation")
  scale <- attr(data_curve, "scale")

  get_procrustes_fit_from_param(data_curve, rot, scale, trans)
}

