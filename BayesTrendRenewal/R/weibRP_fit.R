#' Fit Weibull RP model to data
#'
#' @description Uses Bayesian methods to fit a Weibull renewal process model to a set of volcanoes
#'
#' @param volNames Names of volcanoes to fit (target first)
#' @param VEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param chains Number of MCMC chains to run
#' @param chainLen Length of each chain
#' @param adapt Length of adapt phase for each chain
#'
#' @return Posterior distributions for the parameters of the Weibull RP model (mcmc object)
#' @export
#'
#' @examples ## Fit Weibull model to Taranaki (Target) and analogue volcanoes
#' @examples ## Chain length may need to be longer to get convergence
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples weib.mcmc <- weibRP_fit(volnames, 3, 1000, 3, 5000, 1000)

weibRP_fit <- function(volNames, VEI, period, chains = 3, chainLen, adapt){


  nVol <- length(volNames)
  dat <- readGVP(volNames, minVEI = VEI, period = period, OccurTime = F)

  init_func <- function(){
    return(list(alpha = 1, beta = 1, sigmaA = 0.5, sigmaB = 1))
  }

  m <- jags.model(system.file("/Weib_RP_Mod.txt", package = "BayesTrendRenewal"), data = list(nVol = nVol, x = dat$data[,2] / 1000, n1 = dat$start, n2 = dat$finish, censLimVec = dat$data[,4] / 1000, isCensored = dat$data[,3]),
                  inits = init_func, n.chains = chains, n.adapt = adapt)
  eval(parse(text = paste("out.mcmc <- coda.samples(m, c('alpha', 'beta', 'A', 'B', 'sigmaA', 'sigmaB', 'x[", dat$finish[1], "]'),", chainLen, ")", sep = "")))
  out.mcmc
}
