% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readGVP.R
\name{readGVP}
\alias{readGVP}
\title{Read data from GVP catalog}
\usage{
readGVP(names, minVEI = NULL, period = NULL, OccurTime = T)
}
\arguments{
\item{names}{Names of volcanoes to read (Target volcano first)}

\item{minVEI}{Minimum eruption size (VEI). Default = \code{NULL} e.g. all VEI included}

\item{period}{Time period, e.g. \code{period} = 1000, read data from the last 1000 years. Default = \code{NULL} e.g. all eruptions included}

\item{OccurTime}{If \code{TRUE}, get occurrence times (years since first eruption). If \code{FALSE}, get inter-event times (years). Defaults to \code{TRUE}}
}
\value{
A list containing the data from specified volcanoes
\describe{
\item{\code{$data}}{\itemize{ \item If \code{OccurTime} is \code{TRUE}, a matrix with two columns. First column gives the index for each volcano, second column gives the occurrence times for each volcano in years since first eruption in the specified period.
\item If \code{OccurTime} is \code{FALSE}, a matrix with four columns. First column gives the index for each volcano, second column gives the inter-event times in years or \code{NA} if censored, third column is a binary variable where a 1 indicates a censored inter-event time and a 0 indicates a fully observed inter-event time, fourth column gives 1 + the maximum inter-event time for each volcano for all fully observed inter-event times or the censored time in years for censored inter-event times.}}
\item{\code{$cens}}{A vector containing the censored time for each volcano in years since first eruption in specified period.}
\item{\code{$firsterp}}{A vector containing the first eruption time in specified period for each volcano.}
\item{\code{$start}}{A vector containing the index of the first occurence time or inter-event time for each volcano.}
\item{\code{$finish}}{A vector containing the index of the last occurence time or inter-event time for each volcano.}
\item{\code{$erp}}{A vector contating the number of eruptions in the specified period for each volcano.}
\item{\code{$volcanoName}}{A vector containing the names of each volcano.}
\item{\code{$inter}}{A vector containing the number of inter-event times in specified period for each volcano.}
}
}
\description{
Reads occurrence time data or inter-event time data for eruptions greater than or equal to a specified VEI from one or more volcanoes within a specified time period.
}
\examples{
## Read all occurence time data in the last 10000 years
## from Mt Taranaki and the top potential analogues identified using scheme 1
## Include VEI >= 3 eruptions only
readGVP(c("Taranaki", topLast10000[[1]]), 3, 10000)
}
