#' @title Distance to a smooth curve
#' @description Finds the distance of a discrete open srv curve to a smooth curve
#' @param srv_curve srv transformation of the smooth curve, needs to be vectorised
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @param eps convergence tolerance
#' @import stats
#' @return distance between srv_curve and q

get_distance <- function(srv_curve, s, q, eps = 10*.Machine$double.eps){
  p_integrand <- function(t){sapply(t, function(t) sum(srv_curve(t)^2))}
  p_norm <- integrate(p_integrand,0,1, stop.on.error =FALSe)$value  # should be 1
  q_norm <- sum(q^2 %*% diff(s))  # should be 1
  pq_prod <- sapply(1:(length(s)-1), function(i) {
    integrate(function(t) t(q[,i]) %*% srv_curve(t), s[i], s[i+1], stop.on.error = FALSE)$value
  })
  dist <- sqrt(p_norm + q_norm - 2 * sum(pq_prod))
  dist
}
