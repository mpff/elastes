#' @title Optimal alignment to a smooth curve
#' @description Finds optimal alignment for a discrete open srv curve to a smooth curve
#' @param srv_curve srv transformation of the smooth curve, needs to be vectorised
#' @param s time points for q, first has to be 0, last has to be 1
#' @param q square root velocity vectors, one less than time points in s
#' @param initial_t starting value for the optimisation algorithm
#' @param eps convergence tolerance
#' @import stats
#' @return optimal time points for q, without first value 0 and last value 1,
#' optimal time points have the distance of the observation to the srv_curve as an attribute

find_optimal_t <- function(srv_curve, s, q, initial_t = s, eps = 10*.Machine$double.eps){
  t <- initial_t[2:(length(s)-1)]
  ui <- rbind(diag(1, length(t)), 0) - rbind(0, diag(1, length(t)))
  ci <- c(rep(0, length(t)), -1)

  optim_obj <- constrOptim(t, get_loss, srv_curve = srv_curve, s = s, q = q,
                           ui = ui, ci = ci, grad = get_grad,
                           control = list(reltol = eps))
  t_optim <- c(0, optim_obj$par, 1)
  p_integrand <- function(t){sapply(t, function(t) sum(srv_curve(t)^2))}
  dist <- sqrt(integrate(p_integrand,0,1, stop.on.error = FALSE)$value +
                 sum(q^2%*%diff(s)) + 2*optim_obj$value)
  attr(t_optim, "dist") <- dist
  t_optim
}

get_H <- function(t, srv_curve, s, q){
  # t = time points which need to be optimised
  # srv_curve square root velocity function, needs to be vectorised
  # s = observational time points, first has to be 0, last has to be 1
  # q = observed square root velocity vectors
  t_bounds <- c(0, t, 1)
  h <- function(t,i){
    cross_prod <- q[, i]%*%srv_curve(t)
    (s[i+1] - s[i])*(cross_prod^2)*(cross_prod >= 0)
  }
  H <- sapply(1:ncol(q), function(i){
    integrate(h, i = i, lower = t_bounds[i], upper = t_bounds[i+1],
              rel.tol = 0.1)$value
  })
  attr(H, "t") <- t
  H
}

get_loss <- function(t, srv_curve, s, q){
  H <- get_H(t, srv_curve, s, q)
  assign("H_current", H,envir = parent.frame())
  value <- -sum(sqrt(H))
  return(value)
}

get_grad <- function(t, srv_curve, s, q){
  # t = time points which need to be optimised
  # srv_curve square root velocity function, needs to be vectorised
  # s = observational time points, first has to be 0, last has to be 1
  # q = observed square root velocity vectors
  t_bounds <- c(0, t, 1)
  h <- function(t,i){
    cross_prod <- q[, i]%*%srv_curve(t)
    (s[i+1] - s[i])*(cross_prod^2)*(cross_prod >= 0)
  }
  H <- get0("H_current", envir = parent.frame())
  #check that H actually belongs to t
  if(sum(attributes(H)$t != t) > .Machine$double.eps | is.null(H)){
    H <- sapply(1:ncol(q), function(i){
      integrate(h, i = i, lower = t_bounds[i], upper = t_bounds[i+1],
                rel.tol = 0.1)$value
    })
  }
  h_1 <- sapply(2:ncol(q), function(i) h(t_bounds[i],i-1)/sqrt(H[i-1]))
  h_2 <- sapply(2:ncol(q), function(i) h(t_bounds[i],i)/sqrt(H[i]))
  return(-0.5*(h_1 -h_2))
}


#' Computes elastic distance
#' @noRd
#' @param srv_data_1 srv of curve1
#' @param srv_data_2 srv of curve2
#' @param t_optim optimal parametrisation of curve2

compute_distance<- function(srv_data_1, srv_data_2, t_optim, closed){
  norm_1 <- sum(t(srv_data_1[,-1]^2)%*%diff(c(srv_data_1$t,1)))
  norm_2 <- sum(t(srv_data_2[,-1]^2)%*%diff(c(srv_data_2$t,1)))
  if(closed){
    srv_data_1_extended <- rbind(srv_data_1, srv_data_1, srv_data_1)
    srv_data_1_extended$t <- c(srv_data_1$t - 1, srv_data_1$t, srv_data_1$t + 1)
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1 = srv_data_1_extended, srv_data_2)
  } else {
    cross_prod <- get_loss_discrete(t = t_optim, srv_data_1, srv_data_2)
  }
  sqrt(norm_1 + norm_2 - 2*cross_prod)
}
