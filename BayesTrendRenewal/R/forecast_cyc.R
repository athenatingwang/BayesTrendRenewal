#' Forecast using Cyclic TRP model
#'
#' @description Calculates the posterior distribution for the forecast time for the next eruption at the target volcano using the cyclic TRP model
#'
#' @param out MCMC output for cyclic model
#' @param censorTime Censor time for target volcano (years since first eruption)
#' @param firsterp First eruption time for target volcano (years)
#' @param lasterp Last eruption time for target volcano (years since first eruption)
#' @param Vol Number of volcanoes used in model fitting (target + analogues)
#'
#' @return Posterior distribution for the forecast time (years)
#' @export
#'
#' @examples ## Get median and 95 % credible interval for forecast time at Taranaki using analogue volcanoes
#' @examples ## Chain length may need to be longer to get convergence
#' @examples dat <- readGVP("Taranaki", 3, 1000) # Read taranaki data
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples cyc.mcmc <- cycTRP_fit(volnames, 3, 1000, 3, 5000, 1000)
#' @examples cyc.fore <- forecast_cyc(cyc.mcmc, dat$cens, dat$firsterp, dat$data[dat$erp,2], length(volnames))
#' @examples quantile(cyc.fore, c(0.025, 0.5, 0.975))

forecast_cyc <- function(out, censorTime, firsterp, lasterp, Vol){

  censorTime <- censorTime / 1000
  firsterp <- firsterp / 1000
  lasterp <- lasterp / 1000

  dat <- as.matrix(out)
  nmcmc <- nrow(dat)
  A1mcmc <- dat[,1]
  B1mcmc <- dat[,Vol + 1]
  alphamcmc <- dat[,2*Vol + 1]
  betamcmc <- dat[,2*Vol + 2]
  kappa1mcmc <- dat[,2*Vol + 3]
  mu1mcmc <- dat[,3*Vol + 3]
  omega1mcmc <- dat[,4*Vol + 3]
  phi1mcmc <- dat[,5*Vol + 3]
  rho1mcmc <- dat[,6*Vol + 3]


  ## Calculate shape and rate postior samples for Taranaki
  ratemcmc <- betamcmc * B1mcmc
  shapemcmc <- alphamcmc * A1mcmc

  ## calculate trendfunction parameter for Taranaki
  parammcmc <- kappa1mcmc

  ## P(X < x | X > 2022 - t_n)
  f <- function(x, shape, rate, cens){
    1 - (exp(-(rate*x)^shape) / exp(-(rate*cens)^shape))
  }

  ## Function to solve f(x) = r
  finv <- function(shape, rate, cens, r){

    ftmp <- function(u, r){
      f(u, shape, rate, cens) - r
    }

    uniroot(ftmp, interval = c(0, 100000), r = r)

  }

  forecast <- c()
  for (i in 1:nmcmc){
    r <- runif(1, 0, 1)
    #r = quantile
    cens <- Psi_cyc(censorTime, mu1mcmc[i], phi1mcmc[i], parammcmc[i], omega1mcmc[i], rho1mcmc[i], censorTime) - Psi_cyc(lasterp, mu1mcmc[i], phi1mcmc[i], parammcmc[i], omega1mcmc[i], rho1mcmc[i], censorTime)
    x <- finv(shapemcmc[i], ratemcmc[i], cens, r)$root ## Forecast inter-event (transformed)
    y <- Psi_cyc(lasterp, mu1mcmc[i], phi1mcmc[i], parammcmc[i], omega1mcmc[i], rho1mcmc[i], censorTime) + x ## Add to Time transformed last occurence time
    forecast[i] <- invPsi_cyc(mu1mcmc[i], phi1mcmc[i], parammcmc[i], omega1mcmc[i], rho1mcmc[i], censorTime, y) + firsterp ## Untransform to get actual forecast time

  }

  return(forecast*1000)
}
