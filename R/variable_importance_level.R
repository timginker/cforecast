#' Compute Variable Importance for VAR in Levels
#'
#' Computes the relative importance of each variable in a conditional forecast from a VAR model,
#' taking into account both scale (levels) and variance. Supports both classical VAR models
#' (from the \pkg{vars} package) and Bayesian VAR models (from the \pkg{BVAR} package).
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
#' @importFrom forecast auto.arima
#' @importFrom seasonal seas
#' @importFrom utils head tail
#' @importFrom vars Acoef Bcoef
#' @importFrom methods is
#' @importFrom stats is.ts median resid ts var
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#'
#' @export
#'
variable_importance_level <- function(fit, cond_var, target_var, horizon){


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

  lvls <- variable_weights$`h = 1`

  stds <- NULL

  for (k in 1:ncol(y_hist)) {

    # fit the correct arima model

    arima_model=forecast::auto.arima(y_hist[,k],seasonal = F,
                           max.p = 3,
                           max.q = 3,
                           lambda = 0)


    if(is.numeric(arima_model$model$Delta) && length(arima_model$model$Delta) == 0){

      ndif=0
    } else {
      ndif=1
    }

    arima_model_text=paste0("(",sum(grepl("ar",names(arima_model$coef), ignore.case = TRUE)),
           " ",
           ifelse(ndif==0,0,arima_model$model$Delta[1])," ",
           sum(grepl("ma",names(arima_model$coef), ignore.case = TRUE)),
           ")","(0 0 0)")

    if("drift"%in%names(arima_model$coef) || "intercept" %in%names(arima_model$coef)){

      if(is.ts(y_hist[,k])){


        y_k=y_hist[,k]
      } else{

        y_k=ts(y_hist[,k],frequency = 12,start=2000)

      }

    m_fit <- suppressMessages(seasonal::seas(y_k,
                  regression.aictest = NULL,
                  arima.model = arima_model_text,
                  regression.variables = "const",
                  transform.function = "log",
                  forecast.save = "forecasts",
                  seats.save = "tfd"))
    } else{

      m_fit <- suppressMessages(seasonal::seas(y_k,
                    regression.aictest = NULL,
                    arima.model = arima_model_text,
                    transform.function = "log",
                    forecast.save = "forecasts",
                    seats.save = "tfd"))
    }

    lvls[,k]=as.numeric(head(c(m_fit$series$s12,m_fit$series$tfd),nrow(lvls)))

    stds=c(stds,sqrt(var(resid(m_fit))))


  }




  for (i in seq_len(horizon)) {

    sweights_k <- variable_weights[[i]]

    # Scale each column by the level and standard deviation
    for (j in seq_len(ncol(y_hist))) {
      sweights_k[, j] <- sweights_k[, j] * stds[j]* lvls[,j]
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
