#' State Space Representation
#'
#' Constructs the state space representation of a VAR model. Compatible with both
#' classical VAR models (from the \pkg{vars} package) and Bayesian VAR models
#' (from the \pkg{BVAR} package) for which the state space representation is computed
#' at the median of the parameter distribution.
#'
#' @param fit An object of class \code{varest} (from \pkg{vars}) or \code{bvar} (from \pkg{BVAR})
#' @param p0 diagonal element of the initial state covariance, with default of 1e4
#'
#'
#' @return A list containing the state space system matrices:
#' \itemize{
#'   \item \code{Tt} — transition matrix
#'   \item \code{Zt} — observation matrix
#'   \item \code{GGt} — observation noise covariance
#'   \item \code{HHt} — process noise covariance
#'   \item \code{dt} — deterministic component (e.g., constant)
#'   \item \code{a0} — initial state vector
#'   \item \code{P0} — initial state covariance (diagonal, with default value \code{1e4})
#' }
#'
#' @importFrom vars Acoef Bcoef VAR
#' @importFrom miscTools colMedians
#' @importFrom BVAR companion
#' @importFrom stats median
#'
#' @author Tim Ginker
#' @export
#'
#' @examples
#' library(vars)
#' library(BVAR)
#' data(Canada)
#'
#' # Classical VAR
#' fit_vars <- VAR(Canada, p = 2, type = "const")
#' ss_vars <- state_space_representation(fit_vars)
#'
#' # Bayesian VAR (small example)
#' data <- fred_qd[, c("CPIAUCSL", "UNRATE", "FEDFUNDS")]
#' data <- fred_transform(data, codes = c(5, 5, 1), lag = 4)
#' fit_bvar <- bvar(data, lags = 2, n_draw = 1000L, n_burn = 200L, verbose = FALSE)
#' ss_bvar <- state_space_representation(fit_bvar)

state_space_representation <- function(fit,p0=1e4) {

  if (!inherits(fit, "varest") && !inherits(fit, "bvar")) {
    stop("`fit` must be an object of class 'varest' or 'bvar'.")
  }

  # SS representation for varest objects --------------------

  if (is(fit, "varest")){

  K <- fit$K
  p <- fit$p
  n_state <- K * p

  if (fit$type == "const") {

    # Construct companion matrix
    companion_matrix <- matrix(0, nrow = n_state, ncol = n_state)
    for (j in 0:(p - 1)) {
      companion_matrix[1:K, (j * K + 1):((j + 1) * K)] <- vars::Acoef(fit)[[j + 1]]
    }
    if (p > 1) {
      companion_matrix[(K + 1):n_state, 1:(n_state - K)] <- diag(n_state - K)
    }

    # Augment system for constant
    t1 <- cbind(companion_matrix, rbind(diag(K), matrix(0, ncol = K, nrow = K * (p - 1))))
    t2 <- cbind(matrix(0, nrow = K, ncol = p * K), diag(K))
    Tt <- rbind(t1, t2)

    Zt <- matrix(0, nrow = K, ncol = n_state + K)
    Zt[1:K, 1:K] <- diag(K)

    GGt <- matrix(0, nrow = K, ncol = K)
    ct  <- matrix(0, nrow = K, ncol = 1)

    HHt <- matrix(0, nrow = n_state + K, ncol = n_state + K)
    HHt[1:K, 1:K] <- summary(fit)$covres

    dt <- matrix(0, nrow = n_state + K, ncol = 1)
    a0 <- matrix(0, nrow = n_state + K, ncol = 1)
    P0 <- diag(p0, nrow = n_state + K)

    a0[(p * K + 1):(n_state + K), 1] <- vars::Bcoef(fit)[1:K, K * p + 1]
  }

  if (fit$type == "none") {

    # Construct companion matrix
    companion_matrix <- matrix(0, nrow = n_state, ncol = n_state)
    for (j in 0:(p - 1)) {
      companion_matrix[1:K, (j * K + 1):((j + 1) * K)] <- vars::Acoef(fit)[[j + 1]]
    }
    if (p > 1) {
      companion_matrix[(K + 1):n_state, 1:(n_state - K)] <- diag(n_state - K)
    }

    Tt <- companion_matrix

    Zt <- matrix(0, nrow = K, ncol = n_state)
    Zt[1:K, 1:K] <- diag(K)

    GGt <- matrix(0, nrow = K, ncol = K)
    ct  <- matrix(0, nrow = K, ncol = 1)

    HHt <- matrix(0, nrow = n_state, ncol = n_state)
    HHt[1:K, 1:K] <- summary(fit)$covres

    dt <- matrix(0, nrow = n_state, ncol = 1)
    a0 <- matrix(0, nrow = n_state, ncol = 1)
    P0 <- diag(p0, nrow = n_state)
  }
}
  # SS representation for bvar objects ------------------------


  if (is(fit, "bvar")) {

    K <- dim(fit$beta)[3]
    p <- fit$call$lags
    n_state <- K * p

    companion_matrix <- BVAR::companion(fit,type="quantile")

    # Augment system for constant
    t1 <- cbind(companion_matrix, rbind(diag(K), matrix(0, ncol = K, nrow = K * (p - 1))))
    t2 <- cbind(matrix(0, nrow = K, ncol = p * K), diag(K))
    Tt <- rbind(t1, t2)

    Zt <- matrix(0, nrow = K, ncol = n_state + K)
    Zt[1:K, 1:K] <- diag(K)

    GGt <- matrix(0, nrow = K, ncol = K)
    ct  <- matrix(0, nrow = K, ncol = 1)

    HHt <- matrix(0, nrow = n_state + K, ncol = n_state + K)
    HHt[1:K, 1:K] <- apply(fit$sigma, MARGIN=c(2,3), FUN=median)+diag(10^-6,nrow=ncol(fit$meta$Y))

    dt <- matrix(0, nrow = n_state + K, ncol = 1)
    a0 <- matrix(0, nrow = n_state + K, ncol = 1)
    P0 <- diag(p0, nrow = n_state + K)

    a0[(p * K + 1):(n_state + K), 1] <- miscTools::colMedians(fit$beta[,1,])

  }



  return(list(
    Tt = Tt,
    Zt = Zt,
    GGt = GGt,
    HHt = HHt,
    dt = dt,
    ct = ct,
    a0 = a0,
    P0 = P0
  ))
}
