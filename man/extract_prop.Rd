% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LDA_helper_functions.R
\name{extract_prop}
\alias{extract_prop}
\title{Extract behavior proportion estimates for each track segment}
\usage{
extract_prop(res, ngibbs, nburn, nmaxclust)
}
\arguments{
\item{res}{A list of results returned by \code{\link{cluster_segments}}.
Element \code{theta} stores estimate for behavior proportions for all time
segments.}

\item{ngibbs}{numeric. The total number of iterations of the MCMC chain.}

\item{nburn}{numeric. The length of the burn-in phase.}

\item{nmaxclust}{numeric. A single number indicating the maximum number of
clusters to test.}
}
\value{
A matrix that stores the proportions of each state/cluster (columns)
  per track segment (rows).
}
\description{
Calculates the mean of the posterior for the proportions of each behavior
within track segments. These results can be explored to determine the optimal
number of latent behavioral states.
}
\examples{

\donttest{
#load data
data(tracks.seg)

#select only id, tseg, SL, and TA columns
tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]

#summarize data by track segment
obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))

#cluster data with LDA
res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 1000,
                       nburn = 500, nmaxclust = 7, ndata.types = 2)

#Extract proportions of behaviors per track segment
theta.estim<- extract_prop(res = res, ngibbs = 1000, nburn = 500, nmaxclust = 7)
}


}
