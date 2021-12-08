#' @title Optimal rotation and scaling alignment to a smooth curve
#' @description Finds optimal procrustes alignment for a discrete open srv curve to a smooth curve
#' @param srv_curve.compl srv transformation of the smooth curve, needs to be complex
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @return complex optimal rotation and scaling w for q


find_optimal_w <- function(srv_curve.compl, s, q){
  q.compl <- complex(re = q[,1], im = q[,2])

  # Approximate scalar product
  ints <- sapply(1:(length(s)-1), function(i) {
    re <- integrate(function(t) Re(Conj(q.compl[i]) * srv_curve.compl(t)), s[i], s[i+1], rel.tol = 0.01)$value
    im <- integrate(function(t) Im(Conj(q.compl[i]) * srv_curve.compl(t)), s[i], s[i+1], rel.tol = 0.01)$value
    re + 1i * im
  })
  w_optim <- sum(ints)
  #p_integrand <- function(t){sapply(t, function(t) sum(srv_curve(t)^2))}
  #dist <- sqrt(integrate(p_integrand,0,1, stop.on.error = FALSE)$value +
  #               sum(q^2%*%diff(s)) + 2*optim_obj$value)
  #attr(w_optim, "dist") <- dist
  w_optim
}
