#' Forecast using Gamma or Weibull RP models
#'
#' @description Calculates the posterior distribution for the forecast time for the next eruption at the target volcano using the Gamma or Weibull RP model
#'
#' @param mcmcOut MCMC output for Gamma or Weibull model
#' @param nVol Number of volcanoes used in model fitting (target + analogues)
#' @param lasterp last recorded eruption time at target volcano (years)
#'
#' @return Posterior distribution for the forecast time (years)
#' @export
#'
#' @examples ## Get median and 95 % credible interval for forecast time at Taranaki using analogue volcanoes
#' @examples ## Chain length may need to be longer to get convergence
#' @examples dat <- readGVP("Taranaki", 3, 1000) # Read taranaki data
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples gam.mcmc <- gamRP_fit(volnames, 3, 1000, 3, 5000, 1000)
#' @examples gam.fore <- forecast_RP(gam.mcmc, length(volnames), dat$data[dat$erp,2] + dat$firsterp)
#' @examples quantile(gam.fore, c(0.025, 0.5, 0.975))

forecast_RP <- function(mcmcOut, nVol, lasterp){
  forecast <- as.matrix(mcmcOut[,(nVol * 2) + 5])
  forecast <- (forecast * 1000) + lasterp
  forecast
}
