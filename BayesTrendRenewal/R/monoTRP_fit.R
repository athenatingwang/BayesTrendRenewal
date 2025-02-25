#' Fit monotonic TRP model to data
#'
#' @description Uses Bayesian methods to fit a trend renewal process model with a monotonic trend function to a set of volcanoes
#'
#' @param volNames Names of volcanoes to fit (target first)
#' @param VEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param chains Number of MCMC chains to run
#' @param chainLen Length of each chain
#' @param adapt Length of adapt phase for each chain
#'
#' @return Posterior distributions for the parameters of the monotonic TRP model (mcmc object)
#' @export
#'
#' @examples ## Fit monotonic model to Taranaki (Target) and analogue volcanoes
#' @examples ## Chain length may need to be longer to get convergence
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples mono.mcmc <- monoTRP_fit(volnames, 3, 1000, 3, 5000, 1000)

monoTRP_fit <- function(volNames, VEI, period, chains = 3, chainLen, adapt){


  nVol <- length(volNames)
  dat <- readGVP(volNames, minVEI = VEI, period = period, OccurTime = T)

  init_func <- function(){
    return(list(alpha = 0.5, beta = 0.5, kappa = rep(0.5, nVol), mu = rep(0.5,nVol), sigmaA = 1, sigmaB = 1))
  }

  m <- jags.model(system.file("/Weib_Mod_Trend_Analog.txt", package = "BayesTrendRenewal"), data = list(nVol = nVol, t = dat$data[,2] / 1000, T = dat$cens / 1000, n1 = dat$start, n2 = dat$finish, N = sum(dat$erp - 1)),
                  inits = init_func, n.chains = chains, n.adapt = adapt)
  eval(parse(text = paste("out.mcmc <- coda.samples(m, c('alpha', 'beta', 'kappa', 'mu', 'A', 'B', 'sigmaA', 'sigmaB'),", chainLen, ")", sep = "")))
  out.mcmc
}
