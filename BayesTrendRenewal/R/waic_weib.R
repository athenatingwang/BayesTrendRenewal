#' WAIC Weibull
#'
#' @description Calculates the WAIC value for the Weibull model
#'
#' @param out Output of Weibull model fit
#' @param volnames Names of volcanoes used for model fitting (target first)
#' @param minVEI Minimum eruption size
#' @param period Time period, e.g. period = 1000, use data from the last 1000 years
#' @param target If \code{TRUE}, use only the target volcano data to get WAIC. If False use Analogues also
#'
#' @return WAIC value
#' @export
#'
#' @examples volnames <- c("Taranaki", topLast1000[[1]])
#' @examples weib.mcmc <- weibRP_fit(volnames, 3, 1000, 3, 5000, 1000)
#' @examples weib.waic.analog <- waic_weib(weib.mcmc, volnames, 3, 1000)
#' @examples weib.waic.target <- waic_weib(weib.mcmc, volnames, 3, 1000, target = TRUE)

waic_weib <- function(out, volnames, minVEI, period, target = F){

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
      ppd[k] <- mean(dweibull(inter[j], shape = alpha * A[,i], scale = 1 / (beta * B[,i])))
      pwaic[k] <- var(log(dweibull(inter[j], shape = alpha * A[,i], scale = 1 / (beta * B[,i]))))
      k <- k + 1
    }

    # Censored time
    censored <- cens - occur[length(occur)]

    ppd[k] <- mean(1 - pweibull(censored, shape = alpha * A[,i], scale = 1 / (beta * B[,i])))
    pwaic[k] <- var(pweibull(censored, shape = alpha * A[,i], scale = 1 / (beta * B[,i]), lower.tail = F, log.p = T))
    k <- k + 1

  }

  waic <- -2*sum(log(ppd)) + 2 * sum(pwaic)

  return(waic)
}



