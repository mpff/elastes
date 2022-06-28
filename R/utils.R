#' Calculate the polygon length of a curve
#'
#' @param curve a \code{data.frame} with observed points in each row. Each
#'   variable is one coordinate direction. If there is a variable \code{t},
#'   \code{t_optim} or \code{id}, it is treated as the time parametrization, not
#'   as an additional coordinate.
#' @return The length of \code{curve}, treating it as a polygon.

get_polygon_length <- function(curve) {
  coord_idx <- !(colnames(curve) %in% c("t", "t_optim", "id"))
  curve <- curve[,coord_idx]
  dx <- diff(curve[,1])
  dy <- diff(curve[,2])
  sum(sqrt(dx^2 + dy^2))
}


#' Calculate the center of a curve
#'
#' @param curve a \code{data.frame} with observed points in each row. Each
#'   variable is one coordinate direction. If there is a variable \code{t},
#'   \code{t_optim} or \code{id}, it is treated as the time parametrization, not
#'   as an additional coordinate.
#' @return The average of observed points in \code{curve}.

get_center <- function(curve) {
  coord_idx <- !(colnames(curve) %in% c("t", "t_optim", "id"))
  curve <- curve[,coord_idx]
  colMeans(curve)
}


#' @title Distance to a smooth curve
#' @description Finds the distance of a discrete open srv curve to a smooth curve
#' @param srv_curve srv transformation of the smooth curve, needs to be vectorized
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @param eps convergence tolerance
#' @return distance between srv_curve and q

get_distance <- function(srv_curve, s, q, eps = 10*.Machine$double.eps){
  p_integrand <- function(t){sapply(t, function(t) sum(srv_curve(t)^2))}
  p_norm <- stats::integrate(p_integrand, 0, 1, stop.on.error =FALSE)$value  # should be 1
  q_norm <- sum(q^2 %*% diff(s))  # should be 1
  pq_prod <- sapply(1:(length(s)-1), function(i) {
    stats::integrate(function(t) t(q[,i]) %*% srv_curve(t), s[i], s[i+1], stop.on.error = FALSE)$value
  })
  dist2 <- p_norm + q_norm - 2 * sum(pq_prod)
  if(dist2 < 0) dist2 <- 0
  sqrt(dist2)
}
