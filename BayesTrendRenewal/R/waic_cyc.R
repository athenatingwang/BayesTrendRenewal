#' WAIC Cyclic
#'
#' @description Calculates the WAIC value for the Cyclic TRP model
#'
#' @param out Output of Cyclic model fit
#' @param volnames Names of volcanoes used for model fitting (target first)
#' @param minVEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param target If \code{TRUE}, use only the target volcano data to get WAIC. If \code{FALSE} use Analogues also (Default \code{FALSE})
#'
#' @return WAIC value
#' @export
#'
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples cyc.mcmc <- cycTRP_fit(volnames, 3, 1000, 3, 5000, 1000)
#' @examples cyc.waic.analog <- waic_cyc(cyc.mcmc, volnames, 3, 1000)
#' @examples cyc.waic.target <- waic_cyc(cyc.mcmc, volnames, 3, 1000, target = TRUE)

waic_cyc <- function(out, volnames, minVEI, period, target = F){

  nvol <- length(volnames)
  dat <- readGVP(volnames, minVEI, period)

  A <- as.matrix(out[,1:nvol])
  B <- as.matrix(out[,(nvol + 1) : (nvol*2)])
  alpha <- as.matrix(out[,(nvol * 2 + 1)])
  beta <- as.matrix(out[,(nvol*2 + 2)])
  kappa <- as.matrix(out[,(nvol*2 + 3):(nvol*3 + 2)])
  mu <- as.matrix(out[,(nvol*3 + 3):(nvol*4 + 2)])

  omega <- as.matrix(out[,(nvol*4 + 3):(nvol*5 + 2)])
  phi <- as.matrix(out[,(nvol*5 + 3):(nvol*6 + 2)])
  rho <- as.matrix(out[,(nvol*6 + 3):(nvol*7 + 2)])

  nmcmc <- length(A[,1])

  if (target) {
    ppd <- rep(NA, dat$erp[1])
    pwaic <- rep(NA, dat$erp[1])
  } else {
    ppd <- rep(NA, sum(dat$erp))
    pwaic <- rep(NA, sum(dat$erp))
  }

  k <- 1

  if (target){
    N <- 1
  } else {
    N <- nvol
  }

  for (i in 1:N){

    # Occurence times for vol i
    occur <- dat$data[dat$data[,1] == i,2] / 1000
    # cens time for vol i
    cens <- dat$cens[i] / 1000

    #inter for vol i
    occurT <- matrix(NA, nrow = nmcmc, ncol = length(occur))
    for (j in 1:length(occur)){
      occurT[,j] <- Psi_cyc(occur[j], mu[,i], phi[,i], kappa[,i], omega[,i], rho[,i], cens)
    }
    inter <- matrix(NA, nrow = nmcmc, ncol = length(occur) - 1)
    for (j in 1:nmcmc){
      inter[j,] <- diff(occurT[j,])
    }

    for (j in 1:ncol(inter)){
      ppd[k] <- mean(dweibull(inter[,j], shape = alpha * A[,i], scale = 1 / (beta * B[,i])))
      pwaic[k] <- var(log(dweibull(inter[,j], shape = alpha * A[,i], scale = 1 / (beta * B[,i]))))
      k <- k + 1
    }

    # Censored time
    censored <- Psi_cyc(cens, mu[,i], phi[,i], kappa[,i], omega[,i], rho[,i], cens) - Psi_cyc(occurT[,length(occur)], mu[,i], phi[,i], kappa[,i], omega[,i], rho[,i], cens)

    ppd[k] <- mean(1 - pweibull(censored, shape = alpha * A[,i], scale = 1 / (beta * B[,i])))
    pwaic[k] <- var(pweibull(censored, shape = alpha * A[,i], scale = 1 / (beta * B[,i]), lower.tail = F, log.p = T))
    k <- k + 1

  }

  waic <- -2*sum(log(ppd)) + 2 * sum(pwaic)

  return(waic)
}



