#' WAIC Gamma
#'
#' @description Calculates the WAIC value for the Gamma model
#'
#' @param out Output of Gamma model fit
#' @param volnames Names of volcanoes used for model fitting (target first)
#' @param minVEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param target If TRUE, use only the target volcano data to get WAIC. If False use Analogues also
#'
#' @return WAIC value
#' @export
#'
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples gam.mcmc <- gamRP_fit(volnames, 3, 1000, 3, 5000, 1000)
#' @examples gam.waic.analog <- waic_gam(gam.mcmc, volnames, 3, 1000)
#' @examples gam.waic.target <- waic_gam(gam.mcmc, volnames, 3, 1000, target = TRUE)

waic_gam <- function(out, volnames, minVEI, period, target = F){

  nvol <- length(volnames)
  dat <- readGVP(volnames, minVEI, period)

  A <- as.matrix(out[,1:nvol])
  B <- as.matrix(out[,(nvol + 1) : (nvol*2)])
  alpha <- as.matrix(out[,(nvol * 2 + 1)])
  beta <- as.matrix(out[,(nvol*2 + 2)])

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

    inter <- diff(occur)

    for (j in 1:length(inter)){
      ppd[k] <- mean(dgamma(inter[j], shape = alpha * A[,i], rate = beta * B[,i]))
      pwaic[k] <- var(log(dgamma(inter[j], shape = alpha * A[,i], rate = beta * B[,i])))
      k <- k + 1
    }

    # Censored time
    censored <- cens - occur[length(occur)]

    ppd[k] <- mean(1 - pgamma(censored, shape = alpha * A[,i], rate = beta * B[,i]))
    pwaic[k] <- var(pgamma(censored, shape = alpha * A[,i], rate = beta * B[,i], lower.tail = F, log.p = T))
    k <- k + 1

  }

  waic <- -2*sum(log(ppd)) + 2 * sum(pwaic)

  return(waic)
}



