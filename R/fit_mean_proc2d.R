#' Fitting function for open curves
#' @name fit_mean_proc2d
#' @description Fits an elastic full Procrustes mean for open, planar curves. Is usually called from
#' \code{\link{compute_elastic_proc2d_mean}}.
#' @param srv_data_curves list of \code{data.frame}s with srv vectors in each row.
#' Usually a result of a call to \code{\link{get_srv_from_points}}
#' @param knots set of knots for the mean spline curve
#' @param type if "smooth" linear srv-splines are used which results in a differentiable mean curve
#' if "polygon" the mean will be piecewise linear.
#' @param penalty the penalty to use in the covariance smoothing step. use '-1' for no penalty.
#' @param max_iter maximal number of iterations
#' @param eps the algorithm stops if L2 norm of coefficients changes less
#' @return a \code{list}
#' with entries
#'   \item{type}{"smooth" or "polygon"}
#'   \item{coefs}{\code{coefs} srv spline coefficients of the estimated mean}
#'   \item{knots}{spline knots}
#'   \item{penalty}{penalty usedin the covariance estimation}
#'   \item{t_optims}{optimal parametrisations}
#'   \item{G_optims}{optimal rotations}
#'   \item{b_optims}{optimal scalings}
#'   \item{iter}{number of iterations until convergence}
#'   \item{fit}{a \code{list} containing
#'       \code{gram} the mean basis Gram matrix,
#'       \code{cov_fit} the covariance smoothing objects in the final iteration,
#'       \code{cov_pca} cov coef matrix pca object in the final iteration and
#'       \code{pfit_coefs} the mean basis coefs of smoothed pfits in the final iteration}
#' @importFrom orthogonalsplinebasis SplineBasis GramMatrix
#' @importFrom utils combn
#' @importFrom splines splineDesign

fit_mean_proc2d <- function(srv_data_curves, knots, penalty, max_iter, type, eps){
  # Initial parametrisation, rotation, scaling and coefs.
  t_optims <- lapply(srv_data_curves, function(srv_data_curve){
    c(srv_data_curve$t, 1)
  })
  G_optims <- as.list(rep(0, length(srv_data_curves)))
  b_optims <- as.list(rep(1, length(srv_data_curves)))
  coefs <- 0

  # Get Gram and orthogonal trafo matrix for the mean basis.
  G <- get_gram_matrix(knots, type)
  G.inv <-solve(G)
  chol.G <- chol(G)
  chol.G.inv <- chol(G.inv)

  # Make Integration Grid
  h = 0.01
  arg.grid = seq(0, 1, by=h)

  #iterate procrustes mean fit and warping on procrustes fits
  for (i in 1:max_iter){
    model_data_complex <- get_model_data_complex(t_optims, srv_data_curves, knots, type)

    # build response on s,t-grid per curve
    cov_dat <- lapply(split(model_data_complex, model_data_complex$id), function(x) {
      combs <- combn(1:nrow(x), 2)
      data.frame(
        qq = x$q_m_long[combs[1,]] * Conj(x$q_m_long[combs[2,]]),
        s = x$m_long[combs[1,]],
        t = x$m_long[combs[2,]]
      )
    })
    cov_dat <- do.call(rbind, cov_dat)

    print(dim(cov_dat))

    cov_fit <- smooth_cov_surface(cov_dat, knots, type, penalty)

    beta.mat <- get_complex_coef_matrix(cov_fit, cov_dat)

    # calculate mean from pca on coefficient matrix
    pca <- eigen(chol.G %*% beta.mat %*% t(chol.G))
    coefs.compl <- chol.G.inv %*% Conj(pca$vectors[,1])

    # Normalize mean. ToDo: Mean not normalized correctly for "polygon"
    coefs.norm  <- if(type == "smooth"){
      Re(sqrt(t(Conj(coefs.compl)) %*% G %*% coefs.compl))
    } else {
      mean_eval <- make_design(arg.grid, knots=knots, type = type) %*% coefs.compl
      mm <- Conj(mean_eval) * mean_eval
      sqrt(h * ( sum(mm) - 0.5*mm[1] - 0.5*mm[length(mm)]))
    }
    coefs.compl <- coefs.compl / c(coefs.norm)

    coefs_old <- coefs
    coefs <- as.matrix(data.frame(q_m_long.X1 = Re(coefs.compl), q_m_long.X2 = Im(coefs.compl)))

    # calculate procrustes fits of warped srv data curves, apply rot/scaling to unwarped curves
    srv_procrustes_curves <- lapply(1:length(srv_data_curves), function(j) {

      q <- model_data_complex[model_data_complex$id == j,]

      # Calculate Scalar products qq and qm.
      pfit_prods <- compute_proc2d_alignment(q, coefs.compl, type, knots, h)
      qm <- pfit_prods$qm
      qq <- pfit_prods$qq

      b_optims[[j]] <<- Mod(qm/qq)
      G_optims[[j]] <<- Arg(qm/qq)

      # Calculate procrustes fit of original srv_data_curve
      srv_complex = complex(real = srv_data_curves[[j]][,2], imaginary = srv_data_curves[[j]][,3])
      pfit <- b_optims[[j]] * exp(1.i * G_optims[[j]]) * srv_complex

      # Update srv_data_curve normalization
      #srv_data_curves[[j]] <<- data.frame(t = srv_data_curves[[j]][,1],
      #                                    X1 = 1/b_optims[[j]] * srv_data_curves[[j]][,2])

      # Return SRV Procrustes Fit.
      data.frame(t = srv_data_curves[[j]][,1], X1 = Re(pfit), X2 = Im(pfit))
    })

    #stop if coefficients don't change much anymore
    stop_crit <- sum((coefs - coefs_old)^2)/sum(coefs^2)
    if(stop_crit < eps | max_iter == 0){
      rownames(coefs) <- NULL
      colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
      fit_object <- list(
        "gram" = G,
        "cov_fit" = cov_fit,
        "cov_coef" = beta.mat,
        "cov_pca" = pca,
        "coefs_norm" = coefs.norm,
        "pfit_coefs" = NULL,
        "G_optims" = G_optims,
        "b_optims" = b_optims
      )
      return(list("type" = type, "coefs" = coefs, "knots" = knots,
                  "t_optims" = t_optims, "fit" = fit_object))
    }

    # align parametrization to mean
    t_optims <- get_optimal_t(srv_procrustes_curves, coefs, t_optims, type, knots, eps, i)

  }

  warning("Stopping criteria eps has not been reached! Consider more iterations max_iter")
  rownames(coefs) <- NULL
  colnames(coefs) <- colnames(srv_data_curves[[1]][,-1])
  fit_object <- list(
    "gram" = G,
    "cov_fit" = cov_fit,
    "cov_coef" = beta.mat,
    "cov_pca" = pca,
    "coefs_norm" = coefs.norm,
    "pfit_coefs" = NULL,
    "G_optims" = G_optims,
    "b_optims" = b_optims
  )
  return(list("type" = type, "coefs" = coefs, "knots" = knots,
              "t_optims" = t_optims, "fit" = fit_object))
}

# Get inner and outer knots.
get_knots <- function(knots, type){
  deg <- ifelse(type == "smooth", 1, 0)
  knotl = 1 / ( length(knots) - 1 )  # mean length of a knot
  c(rep(-knotl, deg), knots, rep(1+knotl, deg))
}

# Get b-spline Gram matrix.
get_gram_matrix <- function(knots, type){
  knots = get_knots(knots, type)
  if( type == "polygon"){
    diag(length(knots) - 1)  # Note: ONLY CORRECT FOR EQUIDISTANT KNOTS!!!
  } else {
    osb_smooth = orthogonalsplinebasis::SplineBasis(knots, order = 2)  # degree = order - 1
    orthogonalsplinebasis::GramMatrix(osb_smooth)
  }
}

get_model_data_complex <- function(t_optims, srv_data_curves, knots, type){
  #compute warped srv vectors
  q_m_data <- lapply(1:length(srv_data_curves), function(j){
    old_diff <- diff(c(srv_data_curves[[j]]$t, 1))
    new_diff <- diff(t_optims[[j]])
    as.matrix(srv_data_curves[[j]][,-1])*sqrt(old_diff/new_diff)
  })

  if(type == "polygon"){
    q_m_all <- lapply(1: length(srv_data_curves), function(j){
      q_m_knots <- sapply(knots[-length(knots)], function(knot){
        idx_knot <- findInterval(knot, t_optims[[j]], rightmost.closed = T)
        q_m_data[[j]][idx_knot, ]
      })
      q_knots <- data.frame("t" = knots[-length(knots)], t(q_m_knots))
      q_data <-  data.frame("t" = t_optims[[j]][-length(t_optims[[j]])], q_m_data[[j]])
      data <- rbind(q_knots, q_data)
      unique(data[order(data$t),])
    })
    m <- lapply(q_m_all, function(x){
      c(x$t[-1], 1) - 0.5*diff(c(x$t, 1))
    })
    q_m <- lapply(q_m_all, function(x) x[,-1])
  } else {
    m <- lapply(t_optims, function(t_optim){
      t_optim[-1] - 0.5*diff(t_optim)
    })
    q_m <- q_m_data
  }
  #build datacurve ids
  ids <- lapply(1:length(m), function(j){
    id <- rep(j, length(m[[j]]))
  })
  #convert in long format
  ids_long <- do.call(c, ids)
  m_long <- do.call(c, m)
  q_m_long <- do.call(rbind, q_m)
  q_m_long <- complex(real = q_m_long[,1], imaginary = q_m_long[,2])
  data.frame("id" = ids_long, "m_long" = m_long, "q_m_long" = q_m_long)
}

# creating the design matrix (with outer knots at +- avg knotlegth)
make_design <- function(t, knots, type) {
  deg <- ifelse(type == "smooth", 1, 0)
  knots <- get_knots(knots, type)
  splines::splineDesign(knots = knots, x = t, ord = deg + 1)
}

# Smooth cov surface using mgcv::bam with (skew)-symmetric tensor product p-splines and REML.
smooth_cov_surface <- function(cov_dat, knots, type, penalty){

  # Parameters for covariance smoothing
  cov.m = ifelse(type == "smooth", 0, -1)  # spline degree - 1
  cov.d = penalty # penalty
  cov.fx = ifelse(penalty == -1, TRUE, FALSE)
  cov.knots = get_knots(knots, type)
  cov.k = length(cov.knots) - cov.m - 2  # basis dimension

  # fit covariance surface
  cov_fit <- list()
  cov_fit$re <- mgcv::bam(Re(qq) ~ s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d),
                                     fx = cov.fx, xt = list(skew = FALSE)),
                          data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots))
  cov_fit$im <- mgcv::bam(Im(qq) ~ s(s, t, bs="symm", k = cov.k, m = c(cov.m, cov.d),
                                          fx = cov.fx, xt = list(skew = TRUE)),
                          data = cov_dat, method = "REML", knots=list(s = cov.knots, t = cov.knots))
  cov_fit
}

# get coefficient matrix for pca
get_complex_coef_matrix <- function(cov_fit, cov_dat){
  beta.mat.re <- get_coef_matrix_from_model(cov_fit$re, cov_dat)
  beta.mat.im <- get_coef_matrix_from_model(cov_fit$im, cov_dat)
  matrix(complex(real = as.vector(beta.mat.re), imaginary = as.vector(beta.mat.im)), ncol = ncol(beta.mat.re))
}

get_coef_matrix_from_model <- function(model, cov_dat){
  coefs <- get_unconstrained_basis_coefs(model, cov_dat)
  s_ <- model$smooth[[1]]
  beta <- s_$Z %*% coefs
  F <- s_$bs.dim
  matrix(beta, nrow=F, ncol=F)
}

get_unconstrained_basis_coefs <- function(model, cov_dat){
  cov_dat$qq <- model$y
  s_ <- model$smooth[[1]]
  X <- Predict.matrix(s_, cov_dat)
  print(dim(X))
  X_ <- PredictMat(s_, cov_dat)
  X_ <- cbind(1, X_)
  print(dim(X_))
  D <- solve(crossprod(X), crossprod(X, X_))
  coefs_ <- if(s_$xt$skew == FALSE) {coef(model)} else {c(0, coef(model)[-1])}
  D %*% coefs_
}




