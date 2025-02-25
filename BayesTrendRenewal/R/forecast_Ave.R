#' Model averaged forecasts
#'
#' Get posterior distribution for the model averaged forecasts
#'
#' @param forecasts A Matrix. Each column is the posterior distribution of the forecast time from a different model
#' @param weights Model weights. \code{length(weights)} must be equal to \code{ncol(forecasts)}
#'
#' @return Model averaged posterior distribution for the forecast time
#' @export
#'
#' @examples
#' mod1 <- rgamma(100, 2, 2)
#' mod2 <- rgamma(100, 4, 4)
#' mod3 <- rweibull(100, 2, 1)
#' modAve <- forecast_Ave(cbind(mod1, mod2, mod3), weights = c(0.6, 0.3, 0.1))

forecast_Ave <- function(forecasts, weights){

  out <- rep(NA, nrow(forecasts))

  ind <- sample(1:ncol(forecasts), size = nrow(forecasts), prob = weights, replace = T)

  for (i in 1:nrow(forecasts)){
    out[i] <- forecasts[i,ind[i]]
  }
  out
}

