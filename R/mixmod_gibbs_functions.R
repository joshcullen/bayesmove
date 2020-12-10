#' Internal function to sample latent clusters (for observations)
#'
#' This function samples the latent \emph{z} parameter within the Gibbs sampler.
#' Calls on the \code{rmultinom1} function written in C++. Not for calling
#' directly by users.
#'
#' @param nobs numeric. The total number of rows in the dataset.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param dat A data frame containing only columns of the discretized data
#'   streams for all observations.
#' @param ltheta numeric. A vector of log-transformed estimates for parameter
#'   theta.
#' @param lphi A list containing log-transformed estimates for each data stream
#'   of the phi parameter.
#' @param ndata.types numeric. The number of data streams being analyzed.
#'
#' @return A vector with estimates for \emph{z} for each observation within
#'   \code{dat}.
#'
#'@importFrom stats "runif"
#'
sample.z.mixmod=function(nobs, nmaxclust, dat, ltheta, lphi, ndata.types){
  lprob=matrix(ltheta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA
      tmp=lphi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=0

      #finish calculation
      lprob[,i]=lprob[,i]+tmp
    }
  }
  max1=apply(lprob,1,max)
  lprob=lprob-max1
  tmp=exp(lprob)
  prob=tmp/rowSums(tmp)

  z=rmultinom1(prob=prob, randu=runif(nobs))
  z+1
}
#-----------------------------------

#' Internal function to sample parameter for truncated stick-breaking prior
#'
#' This function samples the latent \emph{v} parameter within the Gibbs sampler.
#' Not for calling directly by users.
#'
#' @param z A vector of latent cluster estimates provided by
#'   \code{\link{sample.z.mixmod}}.
#' @param gamma1 numeric. Hyperparameter for the truncated stick-breaking prior.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#'
#' @return A list with estimates for \emph{v} and \emph{theta} for each of the possible states.
#'
#'
#' @importFrom stats "rbeta"
#'
sample.v.mixmod=function(z, gamma1, nmaxclust){
  tmp=table(z)
  tmp1=rep(0,nmaxclust)
  tmp1[as.numeric(names(tmp))]=tmp
  theta=v=rep(NA,nmaxclust)
  aux=1
  for (i in 1:(nmaxclust-1)){
    seq1=(i+1):nmaxclust
    v[i]=rbeta(1,tmp1[i]+1,sum(tmp1[seq1])+gamma1)
    theta[i]=v[i]*aux
    aux=aux*(1-v[i])
  }
  v[nmaxclust]=1
  theta[nmaxclust]=aux
  list(theta=theta,v=v)
}
#-----------------------------------

#' Internal function to sample bin estimates for each movement variable
#'
#' Estimates values of \emph{phi} matrix for use in characterizing distributions
#' of the movement variables. Not for calling directly by users.
#'
#' @param alpha numeric. A hyperparameter for the Dirichlet distribution.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   data stream. These must be in the same order as the columns within
#'   \code{dat}.
#' @param ndata.types numeric. The number of data streams being analyzed.
#' @param nmat A list based on \code{\link{SummarizeDat}} C++ function to help
#'   with multinomial draws.
#'
#' @return A list of proportion estimates that characterize distributions
#'   (bins) for each data stream and possible behavioral state.
#'
sample.phi.mixmod=function(alpha, nmaxclust, nbins, ndata.types, nmat){
  phi=list()
  for (j in 1:ndata.types){
    tmp2=matrix(NA,nmaxclust,nbins[j])
    for (i in 1:nmaxclust){
      tmp2[i,]=MCMCpack::rdirichlet(1,nmat[[j]][i,]+alpha)
    }
    phi[[j]]=tmp2
  }
  phi
}
#------------------------------------

#' Internal function to calculate the log-likelihood for iteration of mixture model
#'
#' Calculates the log-likelihood of the mixture model based on estimates for
#' \emph{theta} and \emph{phi}.
#'
#' @param phi A list of proportion estimates that characterize distributions
#'   (bins) for each data stream and possible behavioral state.
#' @param theta numeric. A vector of values that sum to one.
#' @param ndata.types numeric. The number of data streams being analyzed.
#' @param dat A data frame containing only columns of the discretized data
#'   streams for all observations.
#' @param nobs numeric. The total number of rows in the dataset.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#'
#' @return A numeric value of the log-likelihood based upon the current values
#'   for \emph{phi} and \emph{theta}.
#'
#'
get.llk.mixmod=function(phi, theta, ndata.types, dat, nobs, nmaxclust){
  prob=matrix(theta,nobs,nmaxclust,byrow=T)
  for (i in 1:nmaxclust){
    for (j in 1:ndata.types){
      #account for NA
      tmp=phi[[j]][i,dat[,j]]
      cond=is.na(dat[,j])
      tmp[cond]=1

      #finish calculation
      prob[,i]=prob[,i]*tmp
    }
  }
  sum(log(rowSums(prob)))
}
#------------------------------------

#' Internal function to sample the gamma hyperparameter
#'
#' @param v numeric. A vector of proportions for each of the possible clusters.
#' @param ngroup numeric. The total number of possible clusters.
#' @param gamma.possib numeric. A vector of possible values that gamma can take
#'   ranging between 0.1 and 1.
#'
#' @return A single numeric value for gamma that falls within
#'   \code{gamma.possib} for calculation of the log-likelihood.
#'
#' @importFrom stats "rmultinom"
#'
sample.gamma.mixmod=function(v, ngroup, gamma.possib){
  #calculate the log probability associated with each possible value of gamma
  cond=v>0.9999999
  v[cond]=0.9999999

  ngamma=length(gamma.possib)
  soma=sum(log(1-v[-ngroup]))
  k=(ngroup-1)*(lgamma(1+gamma.possib)-lgamma(gamma.possib))
  res=k+(gamma.possib-1)*soma

  #exponentiate and normalize probabilities
  res=res-max(res)
  res1=exp(res)
  res2=res1/sum(res1)

  #sample from a categorical distribution
  tmp=rmultinom(1,size=1,prob=res2)
  ind=which(tmp==1)
  gamma.possib[ind]
}
