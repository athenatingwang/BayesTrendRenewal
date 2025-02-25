#' Fit cyclic TRP model to data
#'
#' @description Uses Bayesian methods to fit a trend renewal process model with a cyclic trend function to a set of volcanoes
#'
#' @param volNames Names of volcanoes to fit (target first)
#' @param VEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param chains Number of MCMC chains to run. Default = 3
#' @param chainLen Length of each chain
#' @param adapt Length of adapt phase for each chain
#'
#' @return Posterior distributions for the parameters of the Cyclic TRP model (mcmc object)
#' @export
#'
#' @examples ## Fit cyclic model to Taranaki (Target) and analogue volcanoes
#' @examples ## Chain length may need to be longer to get convergence
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples cyc.mcmc <- cycTRP_fit(volnames, 3, 1000, 3, 5000, 1000)

cycTRP_fit <- function(volNames, VEI, period, chains = 3, chainLen, adapt){


  nVol <- length(volNames)
  dat <- readGVP(volNames, minVEI = VEI, period = period, OccurTime = T)

  init_func <- function(){
    return(list(alpha = 0.5, beta = 0.5, kappa = rep(0.5, nVol), mu = rep(0.5,nVol), phi = rep(0.5, nVol), omega = rep(2, nVol), rho = rep(pi, nVol), sigmaA = 1, sigmaB = 1))
  }

  m <- jags.model(system.file("/Cyc_Mod_Trend_Analog.txt", package = "BayesTrendRenewal"), data = list(nVol = nVol, t = dat$data[,2] / 1000, T = dat$cens / 1000, n1 = dat$start, n2 = dat$finish, N = sum(dat$erp - 1), pi = pi),
                  inits = init_func, n.chains = chains, n.adapt = adapt)
  eval(parse(text = paste("out.mcmc <- coda.samples(m, c('alpha', 'beta', 'kappa', 'mu', 'A', 'B', 'rho', 'omega', 'phi', 'sigmaA', 'sigmaB'),", chainLen, ")", sep = "")))
  out.mcmc
}
