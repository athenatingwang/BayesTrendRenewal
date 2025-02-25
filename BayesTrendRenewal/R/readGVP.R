#' Read data from GVP catalog
#'
#' Reads occurrence time data or inter-event time data for eruptions greater than or equal to a specified VEI from one or more volcanoes within a specified time period.
#'
#'
#' @param names Names of volcanoes to read (Target volcano first)
#' @param minVEI Minimum eruption size (VEI). Default = \code{NULL} e.g. all VEI included
#' @param period Time period, e.g. \code{period} = 1000, read data from the last 1000 years. Default = \code{NULL} e.g. all eruptions included
#' @param OccurTime If \code{TRUE}, get occurrence times (years since first eruption). If \code{FALSE}, get inter-event times (years). Defaults to \code{TRUE}
#'
#' @return
#' A list containing the data from specified volcanoes
#' \describe{
#'    \item{\code{$data}}{\itemize{ \item If \code{OccurTime} is \code{TRUE}, a matrix with two columns. First column gives the index for each volcano, second column gives the occurrence times for each volcano in years since first eruption in the specified period.
#'                                  \item If \code{OccurTime} is \code{FALSE}, a matrix with four columns. First column gives the index for each volcano, second column gives the inter-event times in years or \code{NA} if censored, third column is a binary variable where a 1 indicates a censored inter-event time and a 0 indicates a fully observed inter-event time, fourth column gives 1 + the maximum inter-event time for each volcano for all fully observed inter-event times or the censored time in years for censored inter-event times.}}
#'    \item{\code{$cens}}{A vector containing the censored time for each volcano in years since first eruption in specified period.}
#'    \item{\code{$firsterp}}{A vector containing the first eruption time in specified period for each volcano.}
#'    \item{\code{$start}}{A vector containing the index of the first occurence time or inter-event time for each volcano.}
#'    \item{\code{$finish}}{A vector containing the index of the last occurence time or inter-event time for each volcano.}
#'    \item{\code{$erp}}{A vector contating the number of eruptions in the specified period for each volcano.}
#'    \item{\code{$volcanoName}}{A vector containing the names of each volcano.}
#'    \item{\code{$inter}}{A vector containing the number of inter-event times in specified period for each volcano.}
#'}
#' @export
#'
#' @examples
#' ## Read all occurence time data in the last 10000 years
#' ## from Mt Taranaki and the top potential analogues identified using scheme 1
#' ## Include VEI >= 3 eruptions only
#' readGVP(c("Taranaki", topLast10000[[1]]), 3, 10000)

readGVP <- function(names, minVEI = NULL, period = NULL, OccurTime = T){
  #GVP_Full <- read.csv("GVPdata.csv", skip = 1)
  nVol <- length(names)
  GVP_Analog <- subset(GVP_Full, Volcano.Name %in% names, select = c(2,4,6,9,11,13,15)) ## Select volcanoes of interest and Reduce number of columns
  GVP_Analog <- GVP_Analog[order(GVP_Analog$Volcano.Name),]  ## Order by name

  GVP_Analog$VEI[which(is.na(GVP_Analog$VEI))] = 3    # Set VEI = NA to 3

  if (!is.null(minVEI)) GVP_Analog <- subset(GVP_Analog, VEI >= minVEI) ## remove VEI < 3


  GVP_Analog$Start.Day[which(is.na(GVP_Analog$Start.Day))] = 15  ## Set start day = NA = 15th
  GVP_Analog$Start.Day[which(GVP_Analog$Start.Day == 0)] = 15 ## set day = 0 to 15th
  GVP_Analog$Start.Month[which(is.na(GVP_Analog$Start.Month))] = 7  ## Set start day = NA = 15th
  GVP_Analog$Start.Month[which(GVP_Analog$Start.Month == 0)] = 7 ## set month = 0 to July

  GVP_Analog$erpDate = lubridate::make_date(GVP_Analog$Start.Year, GVP_Analog$Start.Month, GVP_Analog$Start.Day)  ## make column for erruption date

  GVP_Analog$erpDate.dec <- lubridate::decimal_date(GVP_Analog$erpDate)  ## Column for decimal date

  if (!is.null(period)) GVP_Analog <- subset(GVP_Analog, erpDate.dec > 2022 - period)

  GVP_Analog$Volcano.Name  <- factor(GVP_Analog$Volcano.Name, levels = names)
  GVP_Analog <- GVP_Analog[order(GVP_Analog$Volcano.Name), ] # Order data alphabetically with taranaki first

  names = names[which(names %in% GVP_Analog[,1])]
  nVol = length(names)
  ## Numer of eruprions
  n <- c()
  for(i in 1:nVol) n[i] <- length(which(GVP_Analog$Volcano.Name == names[i]))


  nInt <- n - 1

  n2 <- cumsum(n) ## index for last eruption
  n1 <- (n2 - n) + 1 ## index for first eruption

  if (OccurTime == T){

    ## Set up dataframe
    volIdx <- rep(1:nVol, times = n)

    OccurTime <- c()
    firsterp <- c()
    for (i in 1:nVol){
      temp <- GVP_Analog$erpDate.dec[n1[i]:n2[i]]
      temp <- sort(temp)
      firsterp[i] <- min(temp)
      temp <- temp - min(temp)  # set first occur time to 0
      OccurTime <- c(OccurTime, temp)
    }

    dat <- cbind(volIdx, OccurTime)
    censorTime <- 2022 - firsterp
    out <- list(data = dat, cens = censorTime, firsterp = firsterp, start = n1, finish = n2, erp = n, volcanoName = names)
  } else {

    index <- c(); isCensored <- c(); inter <- c(); censLimVec <- c(); lasterp <- c()
    for(i in 1:nVol){
      temp <- GVP_Analog$erpDate.dec[n1[i]:n2[i]]
      temp <- sort(temp)
      lasterp <- append(lasterp, max(temp))
      int <- diff(temp)
      cens <- 2022 - temp[length(temp)]
      inter <- append(inter, c(int, NA))
      isCensored <- append(isCensored, c(rep(0,nInt[i]),1))
      index <- append(index, rep(i, n[i]))
      censLimVec <- append(censLimVec, c(rep(max(int) + 1, nInt[i]), cens))
    }

    dat <- cbind(index, inter, isCensored, censLimVec)
    out <- list(data = dat, start = n1, finish = n2, inter = nInt, volcanoName = names)
  }
  return(out)
}
