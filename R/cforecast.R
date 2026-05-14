#' Conditional Forecast with VAR or BVAR
#'
#' Computes conditional forecasts using the Kalman filter and smoother methods proposed by
#' Clarida and Coyle (1984) and Bańbura et al. (2015). Supports classical VAR models
#' (\code{varest} from the \pkg{vars} package). Future implementations will also include Bayesian VAR models
#' (\code{bvar} from the \pkg{BVAR} package). For now, for BVAR models the function
#' allows analysis at the median of the posterior distribution.
#'
#' @param fit An object of class \code{varest} (from \pkg{vars}) or \code{bvar} (from \pkg{BVAR})
#' @param cond_path A numeric vector or matrix specifying the conditional path for the constrained variables
#' @param cond_var A numeric vector indicating which columns of \code{y} are conditionally constrained
#' @param horizon Optional forecast horizon (number of periods ahead). If \code{NULL}, it is inferred
#'   from the number of rows in \code{cond_path}
#' @param p0 Diagonal element of the initial state covariance matrix. Default is \code{1e4}
#' @param package A character string indicating which backend to use
#'   (\code{"FKF"} or \code{"KFAS"}). Defaults to \code{"FKF"}.
#'   The \code{"KFAS"} backend can be useful when the forecast error
#'   variance matrix is singular or near-singular.
#'
#' @importFrom FKF fkf fks
#' @importFrom KFAS KFS SSModel SSMcustom
#' @importFrom utils tail
#' @importFrom vars Acoef Bcoef
#' @importFrom methods is
#'
#' @references
#'
#' Bańbura, M., Giannone, D., & Lenza, M. (2015). Conditional forecasts and scenario analysis
#'   with vector autoregressions for large cross-sections. \emph{International Journal of Forecasting},
#'   31(3), 739–756.
#'
#' Clarida, R. and D. Coyle (1984). Conditional Projection by Means of Kalman Filtering.
#'   Carnegie-Rochester Conference Series on Public Policy, 20, 247–284.
#'
#' Helske, J. (2017). KFAS: Exponential family state space models in R.
#' \emph{Journal of Statistical Software}, \bold{78}, 1-39.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{forecast} — Conditional forecast matrix (horizon × variables)
#'   \item \code{mse} — Forecast mean squared error (array: K × K × horizon)
#'   \item \code{fkf} — Output of the Kalman smoother (\code{FKF::fks} or \code{KFAS::KFS})
#'   \item \code{ss} — State space representation used in the forecast
#'   \item \code{cond_var} — Indices of constrained variables
#'   \item \code{cond_path} — The conditional path used
#'   \item \code{horizon} — Effective forecast horizon
#'   \item \code{fit} — The fitted VAR/BVAR object
#'   \item \code{y} — Full data matrix (historical + future with conditional constraints)
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
                      package = "FKF") {


  if (!is.null(horizon)) {
    if (!is.numeric(horizon) || length(horizon) != 1 || is.na(horizon) ||
        horizon < 1 || horizon != as.integer(horizon)) {
      stop("`horizon` must be NULL or a positive integer scalar.")
    }
  }


  if (!is.numeric(cond_path)) stop("`cond_path` must be numeric.")

  cond_path <- as.matrix(cond_path)

  if (length(cond_var) == 1 && ncol(cond_path) != 1) {
    # allow vector -> matrix already; but if user passed multi-col accidentally
    stop("With a single `cond_var`, `cond_path` must have exactly 1 column.")
  }

  if (length(cond_var) > 1 && ncol(cond_path) != length(cond_var)) {
    stop("When `cond_var` has length > 1, `cond_path` must have one column per constrained variable.")
  }


  if (!is.null(horizon) && horizon < path_length) {
    stop("`horizon` cannot be smaller than the number of rows in `cond_path`.")
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
    if(package == "FKF"){
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

    }} else if( package == "KFAS"){

      fit_kfas <- SSModel(
        as.matrix(y_all) ~ -1+SSMcustom(
          Z = ss$Zt,
          T = ss$Tt,
          R = diag(dim(ss$Tt)[1]),
          Q = ss$HHt,
          a1 =ss$a0,
          P1 =ss$P0,
          P1inf = matrix(0,nrow = dim(ss$Tt)[1],ncol=dim(ss$Tt)[1])
        ),
        H = ss$GGt
      )

      out <- KFS(fit_kfas, smoothing = c("state", "signal"))

      forecast <- tail(t(ss$Zt%*%t(out$alphahat)),h)


      if(!is.null(colnames(y_hist))){

        colnames(forecast)=colnames(y_hist)
      }

      mse <- array(0,dim=c(dim(forecast)[2],dim(forecast)[2],dim(forecast)[1]))

      for (i in 1:dim(mse)[3]) {

        mse[,,i]=round(t(ss$Zt%*%out$V[,,dim(out$V)[3]-h+i]%*%t(ss$Zt)),10)

      }

      fks_res = out

    }



    return(list(
      forecast = round(forecast,10),
      mse = round(mse,10),
      fkf = fks_res,
      ss = list(a0 = as.numeric(ss$a0), P0 = ss$P0, dt = ss$dt, ct = ss$ct,
                Tt = ss$Tt, Zt = ss$Zt, HHt = ss$HHt, GGt = ss$GGt),
      cond_var = cond_var,
      cond_path = cond_path,
      horizon = h,
      fit = fit,
      y = y_all
    ))




}
