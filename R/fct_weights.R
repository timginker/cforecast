#' Conditional Forecast Variable Weights
#'
#' Computes observation weights used in the Kalman smoother for a given target variable,
#' showing the relative contribution of observed variables to the conditional forecast.
#'
#' @param fit An object of class \code{varest} (from \pkg{vars}) or \code{bvar} (from \pkg{BVAR})
#' @param cond_var A vector indicating which columns of `y` are conditionally constrained
#' @param target_var Column index of the target variable being forecast
#' @param horizon Forecast horizon (number of future periods)
#' @param p0 diagonal element of the initial state covariance, with default of 1e4
#'
#' @importFrom wex wex
#'
#' @return A list of weight matrices for each horizon step
#' @export
#'
#' @examples
#'
#' library(cforecast)
#' library(vars)
#' data(fred_macro)
#' # Fit a  VAR model
#' fit <- VAR(fred_macro[,-1], p = 2, type = "const")
#' # compute observation weights for PCEPILFE given the conditioning on DCOILWTICO
#' fct_weights_fm <- fct_weights(fit = fit, cond_var = 5, target_var = 2, horizon = 1)
#' # note that only DCOILWTICO's weight is nonzero for h=1
#'
fct_weights <- function(fit, cond_var, target_var, horizon, p0 = 1e4) {


  if (!is.numeric(horizon) || length(horizon) != 1 || horizon < 1) {
    stop("`horizon` must be a positive integer.")
  }

  # Prepare historical and placeholder future data
  # varest objects
  if (is(fit, "varest")){
  y_hist <- as.data.frame(fit$y)}
  # bvar objects
  if (is(fit, "bvar")){
    y_hist <- as.data.frame(fit$meta$Y)}

  if (!is.numeric(target_var) || length(target_var) != 1 || target_var < 1 || target_var > ncol(y_hist)) {
    stop("`target_var` must be a valid column index.")
  }


  y_future <- as.data.frame(matrix(0, ncol = ncol(y_hist), nrow = horizon))
  y_future[, -cond_var] <- NA
  colnames(y_future) <- colnames(y_hist)

  y <- rbind(y_hist, y_future)

  # Get state-space representation of the VAR
  ss <- state_space_representation(fit = fit, p0 =p0)

  # Initialize list to store weights
  weights_list <- vector("list", horizon)

  for (k in seq_len(horizon)) {
    idx <- nrow(y_hist) + k

    w_k <- wex::wex(
      Tt  = ss$Tt,
      Zt  = ss$Zt,
      HHt = ss$HHt,
      GGt = ss$GGt,
     # a0  = as.vector(ss$a0),
      P0  = ss$P0,
      yt  = t(y),
      t   = idx
    )

    sweights <- as.data.frame(t(w_k$WtT[target_var, , ]))
    colnames(sweights) <- colnames(y)
    weights_list[[k]] <- sweights
  }

  names(weights_list) <- paste0("h = ", seq_len(horizon))

  return(weights_list)
}
