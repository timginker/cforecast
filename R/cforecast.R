#' Conditional Forecast with VAR
#'
#' Computes conditional forecasts using the methods proposed by Clarida and Coyle (1984) and Banbura et al. (2015).
#'
#' @param fit An object of class 'varest' or 'bvar'
#' @param cond_path A vector or matrix specifying the conditional path
#' @param cond_var A vector indicating which columns of `y` are conditionally constrained
#' @param horizon Forecast horizon. If NULL, the number of rows in `cond_path` is used
#'
#' @import FKF
#' @importFrom utils tail
#' @importFrom vars Acoef Bcoef
#' @importFrom methods is
#'
#' @references Clarida, R. and D. Coyle (1984). Conditional Projection by Means of Kalman Filtering.
#' @references Bańbura, M., Giannone, D. and M. Lenza (2015). Conditional forecasts and scenario analysis with vector autoregressions for large cross-sections. International Journal of Forecasting, 31(3), 739–756.
#'
#' @return A list with the forecast, the forecast MSE, and supporting Kalman state space objects
#' @export
#'
cforecast <- function(fit, cond_path, cond_var, horizon = NULL) {

  if (is(fit, "varest")) {

    K <- fit$K
    p <- fit$p
    n_state <- K * p

    ## Construct companion matrix
    companion_matrix <- matrix(0, nrow = n_state, ncol = n_state)
    for (j in 0:(p - 1)) {
      companion_matrix[1:K, (j * K + 1):((j + 1) * K)] <- vars::Acoef(fit)[[j + 1]]
    }
    if (p > 1) {
      companion_matrix[(K + 1):n_state, 1:(n_state - K)] <- diag(1, nrow = n_state - K)
    }

    ## Kalman system matrices
    Tt <- companion_matrix
    Zt <- matrix(0, nrow = K, ncol = n_state)
    Zt[1:K, 1:K] <- diag(K)
    GGt <- matrix(0, nrow = K, ncol = K)
    ct <- matrix(0, nrow = K, ncol = 1)

    HHt <- matrix(0, nrow = n_state, ncol = n_state)
    HHt[1:K, 1:K] <- summary(fit)$covres

    dt <- matrix(0, nrow = n_state, ncol = 1)
    if (fit$type == "const") {
      dt[1:K] <- vars::Bcoef(fit)[1:K, K * p + 1]
    }

    a0 <- matrix(0, nrow = n_state, ncol = 1)
    P0 <- diag(1e6, nrow = n_state)

    ## Determine horizon and initialize future observations
    cond_path <- as.matrix(cond_path)
    path_length <- nrow(cond_path)
    h <- if (is.null(horizon)) path_length else horizon
    nrows <- max(h, path_length)

    future_obs <- matrix(NA, nrow = nrows, ncol = K)
    for (i in seq_along(cond_var)) {
      future_obs[1:path_length, cond_var[i]] <- cond_path[, i]
    }

    colnames(future_obs) <- colnames(fit$y)
    y_all <- rbind(fit$y, future_obs)

    ## Run Kalman filter and smoother
    fkf_res <- FKF::fkf(
      a0 = as.numeric(a0), P0 = P0, dt = dt, ct = ct,
      Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt, yt = t(y_all)
    )
    fks_res <- FKF::fks(fkf_res)

    ## Extract forecast and MSE
    forecast <- round(tail(t(fks_res$ahatt)[, 1:K], h), 7)
    mse <- fks_res$Vt[1:K, 1:K, (dim(fks_res$Vt)[3] - h + 1):dim(fks_res$Vt)[3]]

    return(list(
      forecast = forecast,
      mse = mse,
      fkf = fks_res,
      ss = list(a0 = as.numeric(a0), P0 = P0, dt = dt, ct = ct,
                Tt = Tt, Zt = Zt, HHt = HHt, GGt = GGt)
    ))
  }

  stop("`fit` must be of class 'varest'. Support for 'bvar' not yet implemented.")
}
