% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixmod_gibbs_functions.R
\name{sample.z.mixmod}
\alias{sample.z.mixmod}
\title{Internal function to sample latent clusters (for observations)}
\usage{
sample.z.mixmod(nobs, nmaxclust, dat, ltheta, lphi, ndata.types)
}
\arguments{
\item{nobs}{numeric. The total number of rows in the dataset.}

\item{nmaxclust}{numeric. A single number indicating the maximum number of
clusters to test.}

\item{dat}{A data frame containing only columns of the discretized data
streams for all observations.}

\item{ltheta}{numeric. A vector of log-transformed estimates for parameter
theta.}

\item{lphi}{A list containing log-transformed estimates for each data stream
of the phi parameter.}

\item{ndata.types}{numeric. The number of data streams being analyzed.}
}
\value{
A vector with estimates for \emph{z} for each observation within
  \code{dat}.
}
\description{
This function samples the latent \emph{z} parameter within the Gibbs sampler.
Calls on the \code{rmultinom1} function written in C++. Not for calling
directly by users.
}
