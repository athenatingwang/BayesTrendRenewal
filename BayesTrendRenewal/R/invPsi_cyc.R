#' Inverse of cyclic trend function
#'
#' @description Uses the cyclic trend function to calculate an occurrence time given a transformed occurrence time
#'
#' @param mu mu parameter
#' @param ph phi parameter
#' @param p kappa parameter
#' @param o omega parameter
#' @param r rho parameter
#' @param c Censored time (time since first eruption)
#' @param b Transformed occurrence time (time since transformed first eruption)
#'
#' @return Occurrence time t (time since first eruption)
#' @export
#'
#' @examples dat <- readGVP("Taranaki", 3, 1000)
#' @examples invPsi_cyc(0.5, 0.5, 4, 3, 3, dat$cens / 1000, 0.5)

invPsi_cyc <- function(mu, ph, p, o, r, c, b){


  f <- function(u, b){
    Psi_cyc(u, mu, ph, p, o, r, c) - b
  }


  uniroot(f, interval = c(0,100000), b = b)$root

}
