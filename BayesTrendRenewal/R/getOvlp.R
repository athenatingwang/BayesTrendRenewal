#' Calculate overlap between posterior distributions
#'
#' @description Compares the posterior distributions of analogue volcanoes to the target volcano by calculating overlap
#'
#' @param x A matrix with n columns. Columns 2:n are compared to column 1
#' @param names A vector of n names, with the first name being the target volcano
#' @param cut_off A value between 0 and 1. Represents the minimum overlap coefficient required for a candidate analogue to be a statistical analogue
#'
#' @return A list. First element is the overlap coefficients between columns 1 and 2:n. Second element is the names of analogues which have an overlap coefficient greater than \code{cut_off}
#' @export
#'
#' @examples ## Get names of analogue volcanoes with A_i overlap > 0.2
#' @examples ## Chain length may need to be longer to get convergence
#' @examples volnames <- c("Taranaki", topLast10000[[1]])
#' @examples gam.mcmc <- gamRP_fit(volnames, 3, 10000, 3, 5000, 1000)
#' @examples Amatrix <- as.matrix(gam.mcmc[,1:length(volnames)])
#' @examples ovlp <- getOvlp(Amatrix, volnames, 0.2)

getOvlp <- function(x, names, cut_off){
  nvol <- ncol(x)
  out <- rep(NA, nvol - 1)

  lim <- range(x)

  den1 <- density(x[,1], from = lim[1] - 1, to = lim[2] + 1)

  area1 <- sfsmisc::integrate.xy(den1$x, den1$y)

  for (i in 2:nvol){

    deni <- density(x[,i], from = lim[1] - 1, to = lim[2] + 1)

    denInt <- pmin(den1$y, deni$y)

    areai <- sfsmisc::integrate.xy(deni$x, deni$y)
    areaInt <- sfsmisc::integrate.xy(den1$x, denInt)

    out[i - 1] <- ((areaInt / area1) + (areaInt / areai)) / 2

  }

  list(ovlp = out, statAnalog = names[2:nvol][out > cut_off])

}
