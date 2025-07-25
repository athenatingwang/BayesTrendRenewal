% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/waic_gam.R
\name{waic_gam}
\alias{waic_gam}
\title{WAIC Gamma}
\usage{
waic_gam(out, volnames, minVEI, period, target = F)
}
\arguments{
\item{out}{Output of Gamma model fit}

\item{volnames}{Names of volcanoes used for model fitting (target first)}

\item{minVEI}{Minimum eruption size}

\item{period}{Time period, e.g. period = 1000, use data from the last 1000 years}

\item{target}{If TRUE, use only the target volcano data to get WAIC. If False use Analogues also}
}
\value{
WAIC value
}
\description{
Calculates the WAIC value for the Gamma model
}
\examples{
volnames <- c("Taranaki", topLast1000[[1]])
gam.mcmc <- gamRP_fit(volnames, 3, 1000, 3, 5000, 1000)
gam.waic.analog <- waic_gam(gam.mcmc, volnames, 3, 1000)
gam.waic.target <- waic_gam(gam.mcmc, volnames, 3, 1000, target = TRUE)
}
