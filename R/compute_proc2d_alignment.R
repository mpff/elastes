#' @title Optimal rotation and scaling alignment to a smooth curve
#' @description Finds optimal rotation and scaling alignment for a discrete open srv curve to a smooth curve
#' @param q complex srv curve with parametrization, needs to be vectorised.
#' The result of a call to \code{get_model_data_complex}
#' @param coefs.compl complex coefficients of smooth curve
#' @param type spline degree
#' @param knots basis knots
#' @return optimal rotation G and scaling b
#' @importFrom stats approx

compute_proc2d_alignment <- function(q, coefs.compl, beta.mat.inv, G, type, knots, h, method = "linear"){

  if (method == "linear"){
    get_intervals <- function(m){
      out <- rep(0, length(m) + 1)
      for(i in 2:length(out)){
        out[i] <- out[i-1] + 2*(m[i-1] - out[i-1])
      }
      return(out)
    }
    m <- get_intervals(q$m_long)

    mean_eval <- make_design(q$m_long, knots = knots, type = type) %*% coefs.compl

    qq <- sum(diff(m) * Conj(q$q_m_long) * q$q_m_long)
    qm <- sum(diff(m) * Conj(q$q_m_long) * mean_eval)
    q_coefs <- NULL

  } else {
    q_B <- make_design(q$m_long, knots, type = type)

    # Estimate coeficients as restricted LSE
    q_coefs <- solve( t(q_B) %*% q_B + Conj(beta.mat.inv)) %*% t(q_B) %*% q$q_m_long

    # Normalize coefs.
    norm <- sqrt(t(Conj(q_coefs)) %*% G %*% q_coefs)
    print(norm)
    q_coefs <- q_coefs / c(norm)

    # Calculate pfit scaling+rotation
    qq <- as.complex(t(Conj(q_coefs)) %*% G %*% q_coefs)
    qm <- as.complex(t(Conj(q_coefs)) %*% G %*% coefs.compl)

  }

  # Calculate G and b
  list("qm" = qm, "qq" = qq, "q_coefs" = q_coefs)
}
