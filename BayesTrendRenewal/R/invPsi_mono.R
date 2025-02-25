#' Inverse of monotonic trend function
#'
#' @description Uses the monotonic trend function to calculate an occurrence time given a transformed occurrence time
#'
#' @param mu mu parameter
#' @param p kappa parameter
#' @param c Censored time (time since first eruption)
#' @param b Transformed occurrence time (time since transformed first eruption)
#'
#' @return Occurrence time t (time since first eruption)
#' @export
#'
#' @examples dat <- readGVP("Taranaki", 3, 1000)
#' @examples invPsi_mono(0.5, 4, dat$cens / 1000, 0.5)

invPsi_mono <- function(mu, p, c, b){

  f <- function(u, b){
    Psi_mono(u, mu, p, c) - b
  }

  uniroot(f, interval = c(0,100000), b = b)$root

}
