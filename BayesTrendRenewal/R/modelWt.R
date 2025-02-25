#' Model Weight
#'
#' Calculates model weights using WAIC values
#'
#' @param waic WAIC values
#'
#' @return Model weights
#' @export
#'
#' @examples
#'   waic <- c(-45, -40, -41, -48)
#'   wt <- modelWt(waic)

modelWt <- function(waic){
  exp(-(waic - min(waic))/2) / sum(exp(-(waic - min(waic))/2))
}
