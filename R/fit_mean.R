#' Mean estimation for open planar curves.
#' @name fit_mean
#' @description Fits an elastic full Procrustes mean for open, planar curves.
#' Is usually called from \code{\link{compute_elastic_shape_mean}}.
#' @param srv_data_curves list of \code{data.frame}s with srv vectors in each row.curves
#' @param knots set of knots for the mean spline curve
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @param penalty the penalty to use in the covariance smoothing step. use '-1' for no penalty.
#' @param pfit_method (experimental) "smooth" or "polygon"
#' @param var_type (experimental) assume "smooth", "constant" or "zero" measurement-error variance along t
#' @param smooth_warp (experimental) controls the weighting of original and smoothed observations
#' over the iterations, if pfit_method == "smooth".
#' @param max_iter maximal number of iterations
#' @param verbose print iterations
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @param cluster a cluster object for use in the \code{bam} call
#' @return a \code{list} with entries
#'   \item{type}{"smooth" or "polygon"}
#'   \item{coefs}{\code{coefs} srv spline coefficients of the estimated mean}
#'   \item{knots}{spline knots}
#'   \item{penalty}{penalty used in the covariance estimation}
#'   \item{distances}{distances to mean}
#'   \item{fit}{a \code{list} containing
#'       \code{t_optims}{optimal parametrizations}
#'       \code{G_optims}{optimal rotations}
#'       \code{b_optims}{optimal scalings}
#'       \code{n_optims}{optimal re-normalization}
#'       \code{n_iter}{number of iterations until convergence}
#'       \code{gram} the mean basis Gram matrix,
#'       \code{cov_fit} the covariance smoothing objects in the final iteration,
#'       \code{cov_pca} cov coef matrix pca object in the final iteration and
#'       \code{pfit_coefs} the mean basis coefs of smoothed pfits in the final iteration}
#' @import sparseFLMM

fit_mean <- function(srv_data_curves, knots, penalty, var_type, pfit_method, max_iter, type, eps, cluster, verbose, smooth_warp){
  # Initial parametrisation, rotation, scaling and coefs.
  t_optims <- lapply(srv_data_curves, function(srv_data_curve){
    c(srv_data_curve$t, 1)
  })
  G_optims <- as.list(rep(0, length(srv_data_curves)))
  b_optims <- as.list(rep(1, length(srv_data_curves)))
  l_optims <- as.list(rep(1, length(srv_data_curves)))
  l_optims_old <- as.list(rep(1, length(srv_data_curves)))
  coefs <- 0
  pfit_coefs <- as.list(rep(0, length(srv_data_curves)))
  distances <- rep(1, length(srv_data_curves))

  # Get Gram and orthogonal trafo matrix for the mean basis.
  G <- get_gram_matrix(knots, type)
  G.inv <-solve(G)
  chol.G <- chol(G)
  chol.G.inv <- chol(G.inv)
  L <- chol(G.inv)
  L.inv <- solve(L)

  # Make Integration Grid
  h = 0.01
  arg.grid = seq(0, 1, by=h)

  #iterate procrustes mean fit and warping on procrustes fits
  for (i in 1:max_iter){

    if(verbose) message(paste0("  Iteration ", i, "..."))

    model_data_complex <- get_model_data_complex(t_optims, srv_data_curves, knots, type)

    # Check for inf in model_data_complex and drop.
    if(any(!is.finite(Re(model_data_complex$q_m_long)))) {
      drops <- model_data_complex[!is.finite(Re(model_data_complex$q_m_long)), ]
      warning(paste("    Warning: Dropping", nrow(drops), "point(s) in mean estimation."))
      model_data_complex <- model_data_complex[is.finite(Re(model_data_complex$q_m_long)), ]
    }

    # Build complex response on s,t-grid per curve
    cov_dat <- lapply(split(model_data_complex, model_data_complex$id), function(x) {
      co <- utils::combn(seq_len(nrow(x)), 2)
      di <- t(matrix(c(seq_len(nrow(x)), seq_len(nrow(x))), ncol=2))  # Include diagonal
      cb <- cbind(di, co)
      combs <- cb[,order(cb[1,], cb[2,])]
      data.frame(
        qq = Conj(x$q_m_long[combs[1,]]) * x$q_m_long[combs[2,]],
        s = x$m_long[combs[1,]],
        t = x$m_long[combs[2,]],
        st = as.numeric(x$m_long[combs[1,]] == x$m_long[combs[2,]])
      )
    })
    cov_dat <- do.call(rbind, cov_dat)

    # Tensor product p-spline smoothing of the complex covariance surface.
    cov_fit <- smooth_cov_surface(cov_dat, knots, type, penalty, var_type, cluster)

    # Get unconstrained basis coef matrix from cov fit.
    beta.mat <- get_complex_coef_matrix(cov_fit)
    beta.mat.inv <- solve(beta.mat)

    # Calculate mean from pca on coefficient matrix
    pca <- eigen(t(L.inv) %*% beta.mat %*% L.inv)
    coefs.compl <- t(L) %*% pca$vectors[,1]

    # Get norm of coefs.
    coefs.norm <- sqrt(t(Conj(coefs.compl)) %*% G %*% coefs.compl)
    coefs.compl <- coefs.compl/c(coefs.norm)

    # Save old coefs, build (x,y) coefs from complex.
    coefs_old <- coefs
    coefs <- as.matrix(data.frame(q_m_long.X1 = Re(coefs.compl), q_m_long.X2 = Im(coefs.compl)))

    # Calculate procrustes fits of the warped srv data curves
    srv_procrustes_curves <- lapply(1:length(srv_data_curves), function(j) {

      q <- model_data_complex[model_data_complex$id == j,]

      # Calculate Scalar products qq and qm (wird noch aufgeräumt!!)
      pfit_prods <- fit_alignment_proc2d(q, type, knots, var_type, coefs.compl, method = pfit_method, cov_fit, pca, L)

      # Get alignment results
      w <- pfit_prods$w  # Conj(z1)
      l <- pfit_prods$l  # length(beta)
      pfit_coefs[[j]] <<- pfit_prods$q_coefs

      # Save scaling and rotation alignment, save re-normalization.
      b_optims[[j]] <<- Mod(w)
      G_optims[[j]] <<- Arg(w)
      l_optims_old[[j]] <<- l_optims[[j]]  # b_optims and l_optims are from diff. iterations!
      l_optims[[j]] <<- l_optims[[j]] * l

      # Apply alignment.
      # SRV curve to complex.
      srv_complex <- complex(real = srv_data_curves[[j]][,2], imaginary = srv_data_curves[[j]][,3])

      # Update normalization of SRV curve.
      srv_complex <- c(1/sqrt(l)) * srv_complex
      srv_data_curves[[j]][,2] <<- Re(srv_complex)
      srv_data_curves[[j]][,3] <<- Im(srv_complex)

      # Blend between smoothed and original curves for warping alignment step.
      if(pfit_method == "smooth"){
        q_coefs_norm <- sqrt(t(Conj(pfit_prods$q_coefs)) %*% G %*% pfit_prods$q_coefs)
        q_coefs_normed <- pfit_prods$q_coefs / c(q_coefs_norm)
        srv_complex_smooth <- make_design(srv_data_curves[[j]][,1], knots, type) %*% q_coefs_normed
        srv_complex <- (1-smooth_warp(i)) * srv_complex + smooth_warp(i) * srv_complex_smooth
      }

      # Calculate rotation aligned and re-normalized SRV curve.
      pfit <- c(w/Mod(w)) * srv_complex  # rotation alignment

      # Return SRV Procrustes Fit.
      srv_procrustes_curve <- data.frame(t = srv_data_curves[[j]][,1], X1 = Re(pfit), X2 = Im(pfit))

      # Update distance to mean
      if(pfit_method == "smooth"){
        distances[[j]] <<- sqrt(1 - Mod(as.matrix(pfit_prods$S)[1,1] + c(w * Conj(w))/l))
      } else {
        distances[[j]] <- get_distance(
          srv_curve = function(t) t(make_design(t, knots=knots, type=type) %*% coefs),
          s = c(srv_procrustes_curve$t, 1),
          q = t(srv_procrustes_curve[,-1]),
          eps = eps*100/i
        )
      }

      srv_procrustes_curve
    })

    #Stop if coefficients didn't change by much.
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if(stop_crit < eps | max_iter == 0){
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      fit_object <- list(
        "n_iter" = i,
        "gram" = G,
        "cov_fit" = cov_fit,
        "cov_coef" = beta.mat,
        "cov_pca" = pca,
        "pfit_coefs" = pfit_coefs,
        "G_optims" = G_optims,
        "b_optims" = b_optims,
        "l_optims" = l_optims_old
      )
      return(list("type" = type, "knots" = knots, "penalty" = penalty,
                  "var_type" = var_type, "pfit_method" = pfit_method, "smooth_warp" = smooth_warp,
                  "coefs" = coefs, "t_optims" = t_optims, "fit" = fit_object, "distances" = distances))
    }

    # align parametrization to mean
    t_optims <- get_optimal_t(srv_procrustes_curves, coefs, t_optims, type, knots, eps, i)

  }

  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
  fit_object <- list(
    "n_iter" = i,
    "gram" = G,
    "cov_fit" = cov_fit,
    "cov_coef" = beta.mat,
    "cov_pca" = pca,
    "pfit_coefs" = pfit_coefs,
    "G_optims" = G_optims,
    "b_optims" = b_optims,
    "l_optims" = l_optims_old
  )
  return(list("type" = type, "knots" = knots, "penalty" = penalty,
              "var_type" = var_type, "pfit_method" = pfit_method, "smooth_warp" = smooth_warp,
              "coefs" = coefs, "t_optims" = t_optims, "fit" = fit_object, "distances" = distances))
}


# Get complex srv data curves, with additional oversampling for identifiability.
get_model_data_complex <- function(t_optims, srv_data_curves, knots, type){

  if (type == "polygon") {
    q_m <- lapply(1:length(srv_data_curves), function(j) {
      curve <- data.frame(t = t_optims[[j]], elasdics::get_points_from_srv(srv_data_curves[[j]]))
      curve_at_knots <- cbind(t = knots, get_evals(curve, t_grid = knots))
      elasdics::get_srv_from_points(curve_at_knots)[, -1, drop = FALSE]
    })
    m <- lapply(srv_data_curves, function(x) {
      knots[-1] - 0.5 * diff(knots)
    })
  }
  else {
    q_m <- lapply(1:length(srv_data_curves), function(j) {
      old_diff <- diff(c(srv_data_curves[[j]]$t, 1))
      new_diff <- diff(t_optims[[j]])
      as.matrix(srv_data_curves[[j]][, -1]) * sqrt(old_diff/new_diff)
    })
    m <- lapply(t_optims, function(t_optim) {
      t_optim[-1] - 0.5 * diff(t_optim)
    })
  }

  #build datacurve ids
  ids <- lapply(1:length(m), function(j){
    id <- rep(j, length(m[[j]]))
  })
  #convert in long format
  ids_long <- do.call(c, ids)
  m_long <- do.call(c, m)
  q_m_long <- do.call(rbind, q_m)

  # Return
  q_m_long <- complex(real = q_m_long[,1], imaginary = q_m_long[,2])
  data.frame("id" = ids_long, "m_long" = m_long, "q_m_long" = q_m_long)
}


# Smooth cov surface using mgcv::bam with (skew)-symmetric tensor product p-splines and REML.
smooth_cov_surface <- function(cov_dat, knots, type, penalty, var_type, cluster){

  # Parameters for covariance smoothing
  cov.m = ifelse(type == "smooth", 0, -1)  # spline degree - 1
  cov.d = penalty # penalty
  cov.fx = ifelse(penalty == -1, TRUE, FALSE)
  cov.knots = get_knots(knots, type)
  cov.k = length(cov.knots) - cov.m - 2  # basis dimension

  # fit covariance surface and measurement-error variance.
  cov_fit <- list()
  if(var_type == "constant"){
    cov_fit$re <- mgcv::bam(
      Re(qq) ~ s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d), fx = cov.fx, xt = list(skew = FALSE)) + st,
      data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots, cluster = cluster)
    )
  } else if(var_type == "smooth") {
    cov_fit$re <- mgcv::bam(
      Re(qq) ~ s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d), fx = cov.fx, xt = list(skew = FALSE))
      + s(t, by = st, bs = "ps", k = cov.k, m = c(cov.m, 1), fx = cov.fx),
      data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots, cluster = cluster)
    )
  } else {
    cov_fit$re <- mgcv::bam(
      Re(qq) ~ s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d), fx = cov.fx, xt = list(skew = FALSE)),
      data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots, cluster = cluster)
    )
  }
  cov_fit$im <- mgcv::bam(
    Im(qq) ~ -1 + s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d), fx = cov.fx, xt = list(skew = TRUE)),
    data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots, cluster = cluster)
  )
  cov_fit
}


# Get coefficient matrix for pca
get_complex_coef_matrix <- function(cov_fit){
  beta.mat.re <- get_coef_matrix_from_model(cov_fit$re)
  beta.mat.im <- get_coef_matrix_from_model(cov_fit$im)
  matrix(complex(real = as.vector(beta.mat.re), imaginary = as.vector(beta.mat.im)), ncol = ncol(beta.mat.re))
}

get_coef_matrix_from_model <- function(model){
  coefs <- get_unconstrained_basis_coefs(model)
  s_ <- model$smooth[[1]]
  beta <- s_$Z %*% coefs
  F <- s_$bs.dim
  matrix(beta, nrow=F, ncol=F)
}

# Calculate the unconstrained basis coefs from constrained mgcv smooth.
get_unconstrained_basis_coefs <- function(model){
  s_ <- model$smooth[[1]]

  # Create toy data that ensures identifiability of trafo mat.
  knots <- s_$knots
  knots_mid <- knots[-length(knots)] + 0.5 * diff(knots)
  knots_mid <- knots_mid[-length(knots_mid)]
  knots_all <- unique(sort(c(knots,knots_mid)))
  toy_dat <- expand.grid(s = knots_all, t = knots_all)

  # Create toy design matrix and calculate trafo matrix.
  X <- mgcv::Predict.matrix(s_, toy_dat)
  X_ <- mgcv::PredictMat(s_, toy_dat)
  X_ <- if(s_$xt$skew == FALSE) {cbind(1, X_)} else {X_}  # ToDo: Check this!
  D <- solve(crossprod(X), crossprod(X, X_))

  # Apply trafo to model coefs.
  coef_idxs <- unique(c(1,s_$first.para:s_$last.para))  # Only intercept and tensor product smooth.
  D %*% stats::coef(model)[coef_idxs]
}

