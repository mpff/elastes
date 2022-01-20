#' @title Optimal rotation and scaling alignment to a smooth curve
#' @description Finds optimal rotation and scaling alignment for a discrete open srv curve to a smooth curve
#' @param q complex srv curve with parametrization, needs to be vectorised.
#' The result of a call to \code{get_model_data_complex}
#' @param coefs.compl complex coefficients of smooth curve
#' @param type spline degree
#' @param knots basis knots
#' @param beta.mat temp
#' @param beta.mat.inv temp
#' @param pen_factor temp
#' @param inv temp
#' @param G tempo
#' @param G.inv temp
#' @param L temp
#' @param pca temp
#' @param method temp
#' @param h temp
#' @return optimal rotation G and scaling b
#' @importFrom stats approx

compute_proc2d_alignment <- function(q, coefs.compl, beta.mat.inv, G, G.inv, L, pca, type, knots, h, method = "linear", pen_factor = 1){

  if (method == "smooth"){
    q_B <- make_design(q$m_long, knots, type = type)
    Lmbd.inv <- diag(1/pca$values)
    V <- pca$vectors

    # Estimate coeficients as restricted LSE
    lambda <- pen_factor
    q_coefs <- t(L) %*% V %*% solve( t(Conj(V)) %*% L %*% t(q_B) %*% q_B %*% t(L) %*% V + lambda * Lmbd.inv ) %*% t(Conj(V)) %*% L %*% t(q_B) %*% q$q_m_long
    #norm <- sqrt(t(Conj(q_coefs)) %*% G %*% q_coefs)
    #q_coefs <- q_coefs/c(norm)

    # Calculate pfit scaling+rotation
    qq <- as.complex(t(Conj(q_coefs)) %*% G %*% q_coefs)
    qm <- as.complex(t(Conj(q_coefs)) %*% G %*% coefs.compl)

  } else {
    m <- get_intervals(q$m_long)

    mean_func <- function(t) {t(make_design(t, knots = knots, type = type) %*% coefs.compl)}
    qm_ints <- sapply(1:(length(m)-1), function(i) {
      re <- integrate(function(s) Re(Conj(q$q_m_long[i]) * mean_func(s)), m[i], m[i+1], rel.tol = 0.01)$value
      im <- integrate(function(s) Im(Conj(q$q_m_long[i]) * mean_func(s)), m[i], m[i+1], rel.tol = 0.01)$value
      re + 1i * im
    })

    qq <- sum(diff(m) * Conj(q$q_m_long) * q$q_m_long)
    qm <- sum(qm_ints)
    q_coefs <- NULL

  }

  # Calculate G and b
  list("qm" = qm, "qq" = qq, "q_coefs" = q_coefs)
}


get_intervals <- function(m){
  out <- rep(0, length(m) + 1)
  for(i in 2:length(out)){
    out[i] <- out[i-1] + 2*(m[i-1] - out[i-1])
  }
  return(out)
}
