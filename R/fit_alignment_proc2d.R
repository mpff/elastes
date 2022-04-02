#' @title Optimal rotation and scaling alignment to a smooth curve
#' @description Finds optimal rotation and scaling alignment for a discrete open srv curve to a smooth curve
#' @param q complex srv curve with parametrization, needs to be vectorised.
#' The result of a call to \code{get_model_data_complex}
#' @param type spline degree
#' @param knots basis knots
#' @param var_type either "smooth" or "constant" measurement error in cov_fit object
#' @param coefs.compl complex coefficients of smooth curve
#' @param method temp
#' @param cov_fit temp
#' @param pca temp
#' @param L temp
#' @return optimal rotation G and scaling b
#' @import stats mgcv

fit_alignment_proc2d <- function(q, type, knots, var_type, coefs.compl, method, cov_fit, pca, L){

  if (method == "smooth") {

    # Get covariance eigen-basis vectors and values.
    V <- pca$vectors
    Lmbd <- pca$values

    # Restrict basis to ensure positive definiteness (see issue #2)
    Lmbd <- Lmbd[pca$values > 0]
    V <- V[ , pca$values >0]

    # Ensure correct Lmbd.inv matrix when only one eigenvalue (see issue #8)
    Lmbd.inv <- if(length(Lmbd) == 1) matrix(1/Lmbd) else diag(1/Lmbd)

    # Construct the measurement error variance
    if(var_type == "zero"){
      T_ <- rep(0, length(q$m_long))
    } else {
      s_index <- ifelse(var_type == "smooth", 2, 1)
      T_ <- predict(cov_fit$re, data.frame(t=q$m_long, s=q$m_long, st=1), type = "terms")[,s_index]
    }

    if(all(T_ > 0)){

      # Construct design matrix in eigen-basis. (Note: L <- chol(G.inv) for Gram-Matrix G of the mean basis)
      q_B <- make_design(q$m_long, knots, type = type)
      E <- q_B %*% t(L) %*% V

      # Construct inverse error variance matrix.
      T.inv <- diag(1/T_)

      # Conditional score covariance.
      S <- solve(t(Conj(E)) %*% T.inv %*% E + Lmbd.inv)

      # Estimate score vector.
      z <- S %*% t(Conj(E)) %*% T.inv %*% q$q_m_long

    } else {

      # Get idxs of pos./neg. part of error variance
      ip <- which(T_ > 0)
      i0 <- which(T_ <= 0)

      if(length(ip) > 0){
        # Construct design matrix in eigen-basis. (Note: L <- chol(G.inv) for Gram-Matrix G of the mean basis)
        q_Bp <- make_design(q$m_long[ip], knots, type = type)
        q_B0 <- make_design(q$m_long[i0], knots, type = type)
        Ep <- q_Bp %*% t(L) %*% V
        E0 <- q_B0 %*% t(L) %*% V

        # Construct inverse error variance matrix (of pos. part).
        Tp.inv <- diag(1/T_[ip])

        # Construct QR-decomp of E0
        QR <- qr(t(E0))
        Q <- qr.Q(QR, complete = T)
        M <- Q[,1:QR$rank]

        # Estimated deterministic part of score vector.
        z0 <- M %*% solve( t(Conj(M)) %*% t(Conj(E0)) %*% E0 %*% M) %*% t(Conj(M)) %*% t(Conj(E0)) %*% q$m_long[i0]

        if(QR$rank < ncol(Q)){
          N <- Q[,(QR$rank + 1):ncol(Q)]

          # Conditional score covariance.
          S <- solve( t(Conj(N)) %*% (t(Conj(Ep)) %*% Tp.inv %*% Ep + Lmbd.inv) %*% N)

          # Estimated random part of score vector.
          zp <- N %*% S %*% t(Conj(N)) %*% ( t(Conj(Ep)) %*% Tp.inv %*% (q$m_long[ip] - Ep %*% z0) - Lmbd.inv %*% z0)
        } else {
          N <- 0; S <- 0; zp <- 0
        }

      } else {
        q_B0 <- make_design(q$m_long[i0], knots, type = type)
        E0 <- q_B0 %*% t(L) %*% V
        QR <- qr(t(E0))
        Q <- qr.Q(QR, complete = T)
        M <- Q[,1:QR$rank]
        z0 <- M %*% solve( t(Conj(M)) %*% t(Conj(E0)) %*% E0 %*% M) %*% t(Conj(M)) %*% t(Conj(E0)) %*% q$m_long[i0]
        N <- 0; S <- 0; zp <- 0
      }

      # Estimate score vector.
      z <- zp + z0

    }

    # Prepare output
    w <- Conj(z[1])
    l <- Re(c(sum(diag(S)) + t(Conj(z)) %*% z))  # length of original curve
    q_coefs <- t(L) %*% V %*% z

  } else {

    m <- get_intervals(q$m_long)

    mean_func <- function(t) {t(make_design(t, knots = knots, type = type) %*% coefs.compl)}
    qm_ints <- sapply(1:(length(m)-1), function(i) {
      re <- integrate(function(s) Re(Conj(q$q_m_long[i]) * mean_func(s)), m[i], m[i+1], rel.tol = 0.01)$value
      im <- integrate(function(s) Im(Conj(q$q_m_long[i]) * mean_func(s)), m[i], m[i+1], rel.tol = 0.01)$value
      re + 1i * im
    })

    w <- sum(qm_ints)
    l <- sum(diff(m) * Conj(q$q_m_long) * q$q_m_long)
    q_coefs <- NULL
    S <- NULL

  }

  # Check length is real and positive. (See issue #2)
  if(Mod(l) != Re(l)){
    stop(paste0("In calculation of curve length: Length is negative or imaginary. l = ", l))
  }

  # Return
  list("w" = w, "l" = Re(l), "q_coefs" = q_coefs, "S" = S)
}


get_intervals <- function(m){
  out <- rep(0, length(m) + 1)
  for(i in 2:length(out)){
    out[i] <- out[i-1] + 2*(m[i-1] - out[i-1])
  }
  return(out)
}
