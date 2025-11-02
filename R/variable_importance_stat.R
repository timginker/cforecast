#' Compute Variable Importance for Stationary VAR
#'
#' Computes the relative importance of each variable in a conditional forecast
#' from a stationary VAR model. Supports both classical VAR models (from the
#' \pkg{vars} package) and Bayesian VAR models (from the \pkg{BVAR} package).
#'
#' @param fit An object of class \code{varest} (from \pkg{vars}) or \code{bvar} (from \pkg{BVAR})
#' @param cond_var A vector indicating which columns of \code{y} are conditionally constrained
#' @param target_var Column index of the target variable being forecasted
#' @param horizon Forecast horizon (number of future periods)
#'
#' @return A list containing two data frames:
#' \itemize{
#'   \item \code{variable_importance} — overall variable importance by horizon
#'   \item \code{marginal_variable_importance} — importance of future observations (marginal)
#' }
#'
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @importFrom utils tail
#' @importFrom stats var
#'
#' @export
variable_importance_stat <- function(fit, cond_var, target_var, horizon) {

  if (!inherits(fit, "varest") && !inherits(fit, "bvar")) {
    stop("`fit` must be an object of class 'varest' or 'bvar'.")
  }

  # varest objects
  if (is(fit, "varest")) {
    y_hist <- as.data.frame(fit$y)
  }
  # bvar objects
  if (is(fit, "bvar")){
    y_hist <- as.data.frame(fit$meta$Y)
  }


  # Compute conditional forecast weights
  variable_weights <- fct_weights(
    fit = fit,
    cond_var = cond_var,
    target_var = target_var,
    horizon = horizon
  )

  variable_importance_list <- vector("list", horizon)
  marginal_importance_list <- vector("list", horizon)

  # Compute standard deviations of each variable in the system
  stds <- sqrt(diag(stats::var(y_hist)))

  for (i in seq_len(horizon)) {

    sweights_k <- variable_weights[[i]]

    # Scale each column by the standard deviation
    for (j in seq_len(ncol(y_hist))) {
      sweights_k[, j] <- sweights_k[, j] * stds[j]
    }

    # Overall Variable Importance
    vi_k <- colSums(abs(sweights_k))
    svi_k <- sum(abs(sweights_k))
    vi_k <- vi_k / svi_k

    variable_importance_list[[i]] <- tibble::tibble(
      horizon  = i,
      variable = names(vi_k),
      share    = vi_k
    )

    # Marginal Variable Importance (future observations only)
    mvi_k <- colSums(abs(tail(sweights_k, horizon)))/ svi_k

    marginal_importance_list[[i]] <- tibble::tibble(
      horizon  = i,
      variable = names(vi_k),
      share    = mvi_k
    )
  }

  variable_importance_df <- dplyr::bind_rows(variable_importance_list)
  marginal_variable_importance_df <- dplyr::bind_rows(marginal_importance_list)

  return(list(
    variable_importance = variable_importance_df,
    marginal_variable_importance = marginal_variable_importance_df
  ))
}
