#' Conditional Forecast with VAR or BVAR
#'
#' Computes conditional forecasts using the Kalman filter and smoother methods proposed by
#' Clarida and Coyle (1984) and Bańbura et al. (2015). Supports classical VAR models
#' (\code{varest} from the \pkg{vars} package) and Bayesian VAR models
#' (\code{bvar} from the \pkg{BVAR} package).
#'
#' @param fit An object of class \code{varest} (from \pkg{vars}) or \code{bvar} (from \pkg{BVAR})
#' @param cond_path A numeric vector or matrix specifying the conditional path for the constrained variables
#' @param cond_var A numeric vector indicating which columns of \code{y} are conditionally constrained
#' @param horizon Optional forecast horizon (number of periods ahead). If \code{NULL}, it is inferred
#'   from the number of rows in \code{cond_path}
#' @param p0 Diagonal element of the initial state covariance matrix. Default is \code{1e4}
#' @param BGL Logical. If \code{TRUE}, implements the method of Banbura et al. (2015). Defaults to \code{FALSE}
#' @param quantiles Optional numeric vector. If BGL is \code{TRUE}, specifies the quantiles of the forecast distribution to compute. By default returns the median.
#'
#' @importFrom FKF fkf fks
#' @importFrom utils tail
#' @importFrom vars Acoef Bcoef
#' @importFrom methods is
#'
#' @references Clarida, R. and D. Coyle (1984). Conditional Projection by Means of Kalman Filtering.
#'   Carnegie-Rochester Conference Series on Public Policy, 20, 247–284.
#'
#' @references Bańbura, M., Giannone, D., & Lenza, M. (2015). Conditional forecasts and scenario analysis
#'   with vector autoregressions for large cross-sections. \emph{International Journal of Forecasting},
#'   31(3), 739–756.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{forecast} — Conditional forecast matrix (horizon × variables)
#'   \item \code{mse} — Forecast mean squared error (array: K × K × horizon)
#'   \item \code{fkf} — Output of the Kalman smoother (\code{FKF::fks})
#'   \item \code{ss} — State space representation used in the forecast
#'   \item \code{cond_var} — Indices of constrained variables
#'   \item \code{cond_path} — The conditional path used
#'   \item \code{horizon} — Effective forecast horizon
#'   \item \code{fit} — The fitted VAR/BVAR object
#'   \item \code{y} — Full data matrix (historical + future with conditional constraints)
#'   \item \code{fct_BGL} — If BGL is \code{TRUE}, returns the forecast quantiles as specified by the \code{quantiles} argument.
#' }
#'
#' @export
#'
#' @examples
#' library(vars)
#' data(Canada)
#' fit <- VAR(Canada, p = 2, type = "const")
#'
#' # Define a conditional path: unemployment at 15 for 3 periods
#' cond_path <- matrix(rep(15, 3), ncol = 1)
#' cond_var <- 4  # 4th column is 'U'
#'
#' conditional_forecast <- cforecast(fit = fit, cond_path = cond_path, cond_var = cond_var)
#' print(conditional_forecast$forecast)

cforecast <- function(fit,
                      cond_path,
                      cond_var,
                      horizon = NULL,
                      p0 = 1e4,
                      BGL = F ,
                      quantiles = NULL) {


  if (!is.null(horizon) && (!is.numeric(horizon) || horizon < 1)) {
    stop("`horizon` must be NULL or a positive integer.")
  }


  if (!inherits(fit, "varest") && !inherits(fit, "bvar")) {
    stop("`fit` must be an object of class 'varest' or 'bvar'.")
  }



  # varest objects
  if (is(fit, "varest")) {
    K=fit$K
    y_hist <- as.data.frame(fit$y)
  }
  # bvar objects
  if (is(fit, "bvar")){
    y_hist <- as.data.frame(fit$meta$Y)
    K=ncol(y_hist)
  }


  if (!is.numeric(cond_var) || any(cond_var < 1 | cond_var > K)) {
    stop("`cond_var` must be valid column indices of the data matrix.")
  }

    ## Create state space representation of the system

    ss <- state_space_representation(fit, p0 = p0)

    ## Determine horizon and initialize future observations
    cond_path <- as.matrix(cond_path)
    path_length <- nrow(cond_path)
    h <- if (is.null(horizon)) path_length else horizon
    nrows <- max(h, path_length)

    future_obs <- matrix(NA, nrow = nrows, ncol = K)
    for (i in seq_along(cond_var)) {
      future_obs[1:path_length, cond_var[i]] <- cond_path[, i]
    }

    colnames(future_obs) <- colnames(y_hist)
    y_all <- rbind(y_hist, future_obs)

    ## Run Kalman filter and smoother
    fkf_res <- FKF::fkf(
      a0 = as.numeric(ss$a0), P0 = ss$P0, dt = ss$dt, ct = ss$ct,
      Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt, yt = t(y_all)
    )
    fks_res <- FKF::fks(fkf_res)

    ## Extract forecast and MSE
    forecast <- tail(t(ss$Zt%*%fks_res$ahatt),h)

    if(!is.null(colnames(y_hist))){

      colnames(forecast)=colnames(y_hist)
    }

    mse <- array(0,dim=c(dim(forecast)[2],dim(forecast)[2],dim(forecast)[1]))

    for (i in 1:dim(mse)[3]) {

      mse[,,i]=round(t(ss$Zt%*%fks_res$Vt[,,dim(fks_res$Vt)[3]-h+i]%*%t(ss$Zt)),10)

    }

    # Implementing the method suggested by Banbura et al. (2015)
    fct_BGL =NULL

    if(BGL){
      K <- fit[["meta"]][["M"]]
      p <- fit[["meta"]][["lags"]]
      n_state <- K * p
      h <- if (is.null(horizon)) path_length else horizon

      fct_arr = array(numeric(0), dim = c(fit$meta$n_save, h, K))

      for (j in 1:fit$meta$n_save) {

        companion_matrix <- matrix(0, nrow = n_state, ncol = n_state)
        companion_matrix[1:K,] <- t(fit$beta[j,-1,])
        companion_matrix[(K + 1):n_state, 1:(n_state - K)] <- diag(n_state - K)

        # Augment system for constant
        t1 <- cbind(companion_matrix, rbind(diag(K), matrix(0, ncol = K, nrow = K * (p - 1))))
        t2 <- cbind(matrix(0, nrow = K, ncol = p * K), diag(K))
        Tt <- rbind(t1, t2)

        Zt <- matrix(0, nrow = K, ncol = n_state + K)
        Zt[1:K, 1:K] <- diag(K)

        GGt <- matrix(0, nrow = K, ncol = K)
        ct  <- matrix(0, nrow = K, ncol = 1)

        HHt <- matrix(0, nrow = n_state + K, ncol = n_state + K)
        HHt[1:K, 1:K] <- fit$sigma[j,,]

        dt <- matrix(0, nrow = n_state + K, ncol = 1)
        a0 <- matrix(0, nrow = n_state + K, ncol = 1)
        P0 <- diag(p0, nrow = n_state + K)

        a0[(p * K + 1):(n_state + K), 1] <- fit$beta[j,1,]

        cond_path = rep(0,5)
        cond_var = 3


        ss = list(
          Tt = Tt,
          Zt = Zt,
          GGt = GGt,
          HHt = HHt,
          dt = dt,
          ct = ct,
          a0 = a0,
          P0 = P0
        )

        ## Determine horizon and initialize future observations
        cond_path <- as.matrix(cond_path)
        path_length <- nrow(cond_path)
        h <- if (is.null(horizon)) path_length else horizon
        nrows <- max(h, path_length)

        future_obs <- matrix(NA, nrow = nrows, ncol = K)
        for (i in seq_along(cond_var)) {
          future_obs[1:path_length, cond_var[i]] <- cond_path[, i]
        }

        colnames(future_obs) <- colnames(y_hist)
        y_all <- rbind(y_hist, future_obs)

        ## Run Kalman filter and smoother
        fkf_res <- FKF::fkf(
          a0 = as.numeric(ss$a0), P0 = ss$P0, dt = ss$dt, ct = ss$ct,
          Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt, yt = t(y_all)
        )
        fks_res <- FKF::fks(fkf_res)

        ## Extract forecast
        forecast <- tail(t(ss$Zt%*%fks_res$ahatt),h)
        fct_arr[j,,]=forecast


      }

      if(is.null(quantiles)){

        fct_BGL <- apply(fct_arr, MARGIN = c(2, 3), FUN = median)
      } else {

        fct_BGL <- apply(fct_arr, c(2, 3), quantile, probs = quantiles)

        fct_BGL <- round(fct_BGL,6)

      }



    }


    return(list(
      forecast = forecast,
      mse = mse,
      fkf = fks_res,
      ss = list(a0 = as.numeric(ss$a0), P0 = ss$P0, dt = ss$dt, ct = ss$ct,
                Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt),
      cond_var = cond_var,
      cond_path = cond_path,
      horizon = h,
      fit = fit,
      y = y_all,
      fct_BGL = fct_BGL
    ))




}
