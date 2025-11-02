#' Compute Conditional Forecast Composition
#'
#' @param x cforecast result
#' @param target_var target variable
#'
#' @returns a data.frame with the decomposition
#' @export
#'
cforecast_composition <- function(x,target_var){

  # Compute forecast variable weights

  forecast_weights <- fct_weights(x$fit,
                                  cond_var = x$cond_var,
                                  target_var = target_var,
                                  horizon = x$horizon)

  # Compute decomposition

  decomposition <- data.frame(matrix(0,nrow = x$horizon,
                                     ncol=ncol(x$y)))

  colnames(decomposition) <- colnames(x$y)

  for (i in seq_len(x$horizon)) {

    decomposition[i,]=colSums(forecast_weights[[i]]*x$y,na.rm = T)

  }

  rownames(decomposition)=names(forecast_weights)

  return(decomposition)

}
