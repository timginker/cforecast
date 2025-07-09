#' Conditional Forecast with VAR
#'
#' Computes conditional forecasts using the methods proposed by Clarida and Coyle (1984) and Banbura et al. (2015).
#'
#' @param fit An object of class 'varest' or bvar
#' @param cond_path A vector or a matrix with conditional path
#' @param cond_var A vector with constrained columns in y
#' @param horizon forecast horizon
#'
#' @import FKF
#' @importFrom utils tail
#' @importFrom vars Acoef Bcoef
#' @importFrom methods is
#'
#' @references Clarida, R. and D. Coyle (1984). Conditional Projection by Means of Kalman Filtering.
#' @references Bańbura, M., Giannone, D. and M. Lenza (2015) . Conditional forecasts and scenario analysis with vector autoregressions for large cross-sections. International Journal of forecasting, 31(3), pp.739-756.
#'
#' @returns A list with forecasts and the associated MSE
#' @export
#'
cforecast<-function(fit,
                    cond_path,
                    cond_var,
                    horizon=NULL){

  # if(!(class(fit)%in%c('varest','bvar'))){
  #   stop("fit argument should be of class varest or bvar")
  # }

  if(is(fit,"varest")){

    # creating companion matrix

    companion=matrix(0,nrow = fit$p*fit$K,ncol=fit$p*fit$K)

    for(j in 0:(fit$p-1)){

      companion[1:fit$K,(j*fit$K+1):((j+1)*fit$K)]=vars::Acoef(fit)[[j+1]]

    }

    if(fit$p>1){

      companion[(fit$K+1):(fit$K*fit$p),1:(fit$K*(fit$p-1))]=diag(1,nrow = (fit$K*(fit$p-1)))

    }

    Tt=companion

    Zt=matrix(0,nrow = fit$K, ncol= fit$K*fit$p)

    Zt[1:fit$K,1:fit$K]=diag(1,nrow = fit$K)

    GGt=matrix(0,nrow = fit$K,ncol = fit$K)

    ct=matrix(0,nrow = fit$K,ncol=1)

    HHt=matrix(0,nrow = fit$p*fit$K,ncol = fit$p*fit$K)

    HHt[1:fit$K,1:fit$K]=summary(fit)$covres

    dt=matrix(0,nrow = fit$p*fit$K,ncol=1)

    if(fit$type=="const"){

      dt[1:fit$K]=vars::Bcoef(fit)[1:fit$K,(fit$K*fit$p+1)]

    }

    a0=matrix(0,nrow = fit$p*fit$K, ncol=1)

    P0=diag(10^6,nrow = fit$p*fit$K)

    # adding rows

    if(is.null(horizon)){

      horizon=nrow(as.matrix(cond_path))

      future_obs=matrix(NA,nrow = horizon,ncol=fit$K)


    }else{

      future_obs=matrix(NA,nrow = max(nrow(as.matrix(cond_path)),horizon),ncol=fit$K)

    }

    for (i in 1:length(cond_var)) {

      future_obs[1:nrow(as.matrix(cond_path)),cond_var[i]]=as.matrix(cond_path)[,i]

    }

    colnames(future_obs)=colnames(fit$y)

    yfin=rbind(fit$y,future_obs)


    res=FKF::fkf(a0 = as.numeric(a0),
                 P0 =  P0,
                 dt =dt,
                 ct = ct,
                 Tt = Tt,
                 Zt = Zt,
                 HHt = HHt,
                 GGt = GGt,
                 yt = t(yfin))


    fct=round(utils::tail(t(res$att)[,1:fit$K],horizon),7)

    mse=res$Ptt[1:fit$K,1:fit$K,(dim(res$Ptt)[3]-horizon+1):dim(res$Ptt)[3]]


  }


  return(list(forecast=fct,mse=mse,fkf=res))


}
