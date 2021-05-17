#' @title Finds optimal alignment for discrete open curves
#' @description Finds optimal aligned time points for srv curve q to  srv curve p using
#' coordinate wise optimisation.
#' @param srv_procrustes_curves scaling and rotation aligned srv curves
#' @param coefs mean coefficients
#' @param t_optims current optimal parametrisation
#' @param type "smooth" or "polygon"
#' @param knots mean basis knots
#' @param eps convergence tolerance
#' @param i current iteration
#' @return optimal time points for srv_data_curves, without first value 0 and last value 1
#' optimal time points have the distance of the observation to the srv_curve as an attribute

get_optimal_t <- function(srv_procrustes_curves, coefs, t_optims, type, knots, eps, i){
  if(type == "polygon" ){
    t_optims <- lapply(1:length(srv_procrustes_curves), function(j){
      t_optim <- find_optimal_t_discrete(r = knots,
                                         p = t(coefs),
                                         s = c(srv_procrustes_curves[[j]]$t, 1),
                                         q = t(srv_procrustes_curves[[j]][,-1]),
                                         initial_t = t_optims[[j]],
                                         eps = eps*100/i)
      attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
      attr(t_optim, "dist") <- NULL
      t_optim
    })
  } else {
    pfun <- function(t){
      t(make_design(t, knots = knots, type = type) %*% coefs)
    }
    t_optims <- lapply(1:length(srv_procrustes_curves), function(j){
      t_optim <- find_optimal_t(srv_curve = pfun,
                                s = c(srv_procrustes_curves[[j]]$t, 1),
                                q = t(srv_procrustes_curves[[j]][,-1]),
                                initial_t = t_optims[[j]],
                                eps = eps*100/i)
      attr(t_optim, "dist_to_mean") <- attr(t_optim, "dist")
      attr(t_optim, "dist") <- NULL
      t_optim
    })
  }
  t_optims
}
