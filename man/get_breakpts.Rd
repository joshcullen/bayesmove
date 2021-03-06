% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation_helper_functions.R
\name{get_breakpts}
\alias{get_breakpts}
\title{Extract breakpoints for each animal ID}
\usage{
get_breakpts(dat, MAP.est)
}
\arguments{
\item{dat}{A list of lists where animal IDs are separated as well as the
breakpoints estimated for each iteration of the MCMC chain. This is stored
within \code{breakpts} of model results returned after running
\code{\link{segment_behavior}}.}

\item{MAP.est}{numeric. A vector of values at which the maximum a posteriori
(MAP) estimate was identified for each of the animal IDs as returned by
\code{\link{get_MAP}}. These must be in the same order as the data for the
IDs supplied to \code{segment_behavior()}.}
}
\value{
A data frame where breakpoints are returned per animal ID within each
  row. For animal IDs that have fewer breakpoints than the maximum number
  that were estimated, \code{NA} values are used as place holders for these
  breakpoints that do not exist.
}
\description{
Extract breakpoints for each animal ID
}
\examples{
\donttest{
#load data
data(tracks.list)

#subset only first track
tracks.list<- tracks.list[1]

#only retain id and discretized step length (SL) and turning angle (TA) columns
tracks.list2<- purrr::map(tracks.list,
                   subset,
                  select = c(id, SL, TA))


set.seed(1)

# Define model params
alpha<- 1
ngibbs<- 1000
nbins<- c(5,8)

#future::plan(future::multisession)  #run all MCMC chains in parallel
dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
                           alpha = alpha)


# Determine MAP iteration for selecting breakpoints and store breakpoints
MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)
}

}
