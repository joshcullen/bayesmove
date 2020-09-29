#' Internal function to sample latent clusters
#'
#' This function samples the latent \emph{z} parameter within the Gibbs sampler.
#' Calls on the \code{SampleZAgg} function written in C++. Not for calling
#' directly by users.
#'
#' @param ntsegm numeric. The total number of time segments from all animal IDs.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{y}.
#' @param y A list where each element stores separate aggregated count data per
#'   bin per time segment for each movement variable being analyzed. These are
#'   stored as matrices.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param phi A list where each element stores separate proportions per bin per
#'   time segment for each movement variable.
#' @param ltheta A matrix storing the log-transformed values from the
#'   \code{theta} parameter.
#' @param zeroes A list of arrays that contain only zero values which are three
#'   dimensional (ntsegm,nbins[i],nmaxclust).
#' @param ndata.types numeric. The number of data streams being analyzed.
#'
#' @return A list with estimates for \emph{z} where the number of elements is equal to
#'   the number of movement variables.
#'
#'
#'
#'
sample.z=function(ntsegm, nbins, y, nmaxclust, phi, ltheta, zeroes, ndata.types){
  z.agg=list()
  for (i in 1:ndata.types){
    tmp=SampleZAgg(ntsegm=ntsegm,b1=nbins[i],y1=y[[i]],nmaxclust=nmaxclust,
                   lphi1=log(phi[[i]]), ltheta=ltheta,zeroes=zeroes[[i]])
    z.agg[[i]]=tmp$Z1Agg
  }
  z.agg
}
#-----------------------------------

#' Internal function to sample parameter for truncated stick-breaking prior
#'
#' This function samples the latent \emph{v} parameter within the Gibbs sampler.
#' Calls on the \code{CumSumInv} function written in C++. Not for calling
#' directly by users.
#'
#' @param z.agg A list of latent cluster estimates provided by
#'   \code{\link{sample.z}}.
#' @param gamma1 numeric. Hyperparameter for the truncated stick-breaking prior.
#' @param ntsegm numeric. The total number of time segments from all animal IDs.
#' @param ndata.types numeric. The number of data streams being analyzed.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#'
#' @return A matrix with estimates for \emph{v} for each of the number of time
#'   segments and possible states.
#'
#'
#'
#' @importFrom stats "rbeta"
#'
sample.v=function(z.agg, gamma1, ntsegm, ndata.types, nmaxclust){
  soma.fim=matrix(0,ntsegm,nmaxclust)
  cumsum.fim=matrix(0,ntsegm,nmaxclust-1)
  for (i in 1:ndata.types){
    soma=apply(z.agg[[i]],c(1,3),sum)
    soma.fim=soma.fim+soma
    tmp=CumSumInv(ntsegm=ntsegm,nmaxclust=nmaxclust,z=soma)
    cumsum.fim=cumsum.fim+tmp[,-1]
  }

  v=matrix(NA,ntsegm,nmaxclust-1)
  for (i in 1:(nmaxclust-1)){
    v[,i]=rbeta(ntsegm,soma.fim[,i]+1,cumsum.fim[,i]+gamma1)
  }
  cbind(v,1)
}
#-----------------------------------

#' Internal function to calculate theta parameter
#'
#' Calculates values of \emph{theta} matrix within Gibbs sampler. Not for calling
#' directly by users.
#'
#' @param v A matrix returned by \code{\link{sample.v}}
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param ntsegm numeric. The total number of time segments from all animal IDs.
#'
#' @return A matrix of proportion estimates that represent proportions of
#'   different behavioral states per time segment.
#'
#'
#'
#'
#'
get.theta=function(v, nmaxclust, ntsegm){
  theta=matrix(NA,ntsegm,nmaxclust)
  theta[,1]=v[,1]
  prod=rep(1,ntsegm)
  for (i in 2:nmaxclust){
    prod=prod*(1-v[,i-1])
    theta[,i]=v[,i]*prod
  }
  theta
}
#-----------------------------------

#' Internal function to sample bin estimates for each movement variable
#'
#' Estimates values of \emph{phi} matrix for use in characterizing distributions
#' of the movement variables. Not for calling directly by users.
#'
#' @param z.agg A list of latent cluster estimates provided by
#'   \code{\link{sample.z}}.
#' @param alpha numeric. A hyperparameter for the Dirichlet distribution.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{y}.
#' @param ndata.types numeric. The number of data streams being analyzed.
#'
#' @return A matrix of proportion estimates that characterize distributions
#'   (bins) for each movement variable and possible behavioral state.
#'
#'
#'
#'
#'
sample.phi=function(z.agg, alpha, nmaxclust, nbins, ndata.types){
  phi=list()
  for (j in 1:ndata.types){
    soma=t(apply(z.agg[[j]],2:3,sum))
    tmp=matrix(NA,nmaxclust,nbins[j])
    for (i in 1:nmaxclust){
      tmp[i,]=MCMCpack::rdirichlet(1,soma[i,]+alpha)
    }
    phi[[j]]=tmp
  }
  phi
}

