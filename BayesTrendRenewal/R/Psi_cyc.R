#' Cyclic trend function
#'
#' @description Function used to transform occurrence times using the cyclic trend function
#'
#' @param t Occurrence time (time since first eruption)
#' @param mu mu parameter
#' @param phi phi parameter
#' @param kappa kappa parameter
#' @param omega omega parameter
#' @param rho rho parameter
#' @param tt Censored time (time since first eruption)
#'
#' @return Transformed occurrence time (time since transformed first eruption)
#' @export
#'
#' @examples dat <- readGVP("Taranaki", 3, 1000)
#' @examples Psi_cyc(dat$data[,2] / 1000, mu = 0.5, phi = 0.5, kappa = 4, omega = 3, rho = 3 ,dat$cens / 1000)

Psi_cyc <- function(t, mu, phi, kappa, omega, rho, tt){
  return((mu * t) + (1 - mu) * ((phi*tt*(t/tt)^kappa) + (1 - phi)*(t + (tt/(omega*pi)) * sin(omega*pi*(t/tt) + rho) * sin(omega*pi*(t/tt)))))
}
