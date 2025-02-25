#' Monotonic trend function
#'
#' @description Function used to transform occurrence times using the monotonic trend function
#'
#' @param t Occurrence time (time since first eruption)
#' @param mu mu parameter
#' @param kappa kappa parameter
#' @param tt Censored time (time since first eruption)
#'
#' @return Transformed occurrence time (time since transformed first eruption)
#' @export
#'
#' @examples dat <- readGVP("Taranaki", 3, 1000)
#' @examples Psi_mono(dat$data[,2] / 1000, mu = 0.5, kappa = 4, dat$cens / 1000)

Psi_mono <- function(t, mu, kappa, tt){
  return(mu * t + (1 - mu) * tt * (t/tt)^kappa)
}
