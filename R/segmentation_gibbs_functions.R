#' Internal function that calculates the sufficient statistics for the
#' segmentation model
#'
#' An internal function that calculates the sufficient statistics to be used
#' within the reversible-jump MCMC Gibbs sampler called by
#' \code{link{samp_move}}.
#'
#' @param breakpt numeric. A vector of breakpoints.
#' @param dat A matrix that only contains columns storing discretized data for
#'   each of the movement variables.
#' @param max.time numeric. The number of of the last observation of \code{dat}.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{dat}.
#' @param ndata.types numeric. The length of \code{nbins}.
#'
#' @return Returns the sufficient statistics associated with the provided
#'   breakpoints for a given animal ID.
#'
#'
#'
#'
#'
get_summary_stats=function(breakpt, dat, max.time, nbins, ndata.types){

  breakpt1<- c(1, breakpt, max.time)
  n<- length(breakpt1)

  #get summarized results
  res<- list()
  for (j in 1:ndata.types){
    #initialize matrix
    res[[j]]<- matrix(0, n-1, nbins[j])

    #get summary for each interval
    for (i in 2:n){
      if (i < n)
        ind<- breakpt1[i-1]:(breakpt1[i]-1)
      if (i == n)
        ind<- breakpt1[i-1]:(breakpt1[i])

      tmp<- dat[ind,j]
      tmp1<- table(tmp)
      ind<- as.numeric(names(tmp1))
      res[[j]][i-1,ind]<- tmp1
    }
  }

  res
}
#---------------------------------------------

#' Internal function that calculates the log marginal likelihood of each model
#' being compared
#'
#' An internal function that is used to calculate the log marginal likelihood of
#' models for the current and proposed sets of breakpoints. Called within
#' \code{\link{samp_move}}.
#'
#' @param alpha numeric. A single value used to specify the hyperparameter for
#'   the prior distribution. A standard value for \code{alpha} is typically 1,
#'   which corresponds with a vague prior on the Dirichlet distribution.
#' @param summary.stats A matrix of sufficient statistics returned from
#'   \code{\link{get_summary_stats}}.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable.
#' @param ndata.types numeric. The length of \code{nbins}.
#'
#' @return The log marginal likelihood is calculated for a model with a given
#'   set of breakpoints and the discretized data.
#'
#'
#'
#'
#'
log_marg_likel=function(alpha, summary.stats, nbins, ndata.types){

  #get ratios
  probs<- rep(NA, ndata.types)
  for (i in 1:ndata.types){
    lnum<- rowSums(lgamma(alpha + summary.stats[[i]]))
    lden<- lgamma(nbins[i] * alpha+rowSums(summary.stats[[i]]))
    p2<- sum(lnum) - sum(lden)
    p1<- nrow(summary.stats[[i]]) * (lgamma(nbins[i] * alpha) - nbins[i] * lgamma(alpha))
    probs[i]<- p1 + p2
  }

  sum(probs)
}
#---------------------------------------------

#' Internal function for the Gibbs sampler within the reversible-jump MCMC
#' algorithm
#'
#' This is RJMCMC algorithm that drives the proposal and selection of
#' breakpoints for the data based on the difference in log marginal likelihood.
#' This function is called within \code{\link{behav_gibbs_sampler}}.
#'
#' @param breakpt numeric. A vector of breakpoints.
#' @param max.time numeric. The number of of the last observation of \code{dat}.
#' @param dat A matrix that only contains columns storing discretized data for
#'   each of the movement variables used within \code{\link{get_summary_stats}}.
#' @param alpha numeric. A single value used to specify the hyperparameter for
#'   the prior distribution. A standard value for \code{alpha} is typically 1,
#'   which corresponds with a vague prior on the Dirichlet distribution.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{dat}.
#' @param ndata.types numeric. The length of \code{nbins}.
#'
#' @return The breakpoints and log marginal likelihood are retained from the
#'   selected model from the Gibbs sampler and returned as elements of a list.
#'   This is performed for each iteration of the MCMC algorithm.
#'
#'
#'
#' @importFrom stats "runif"
#'
samp_move=function(breakpt, max.time, dat, alpha, nbins, ndata.types){

  breakpt.old<- breakpt
  p<- length(breakpt)
  new.brk<- sample(2:max.time, size=1)  #don't propose a new.brk at brkpt = 1
  brk.augmented<- sort(unique(c(breakpt.old, new.brk)))

  p0<- 1
  rand1<- runif(1)
  if (p == 1) {
    #birth
    if (rand1 < 1/2){
      breakpt.new<- brk.augmented
      p0<- 2/3 #death prob 2 -> 1 is (1/3) and birth prob 1 -> 2 is 1/2.
    }
    #swap
    if (rand1 > 1/2) breakpt.new=new.brk
  }
  if (p > 1) {
    #birth
    if (rand1 < 1/3) {
      breakpt.new<- brk.augmented
    }
    #death
    if (rand1 > 1/3 & rand1 < 2/3) {
      ind<- sample(1:length(breakpt.old), size=1)
      breakpt.new<- breakpt.old[-ind]
      if (p==2) p0<- 3/2 #birth prob from 1 -> 2 is 1/2 and death prob from 2 -> 1 is 1/3
    }
    #swap
    if (rand1 > 2/3) {
      ind<- sample(1:length(breakpt.old), size=1)
      breakpt.new<- sort(unique(c(breakpt.old[-ind], new.brk)))
    }
  }

  #get sufficient statistics
  stats.old<- get_summary_stats(breakpt = breakpt.old, max.time = max.time, dat = dat,
                              nbins = nbins, ndata.types = ndata.types)
  stats.new<- get_summary_stats(breakpt = breakpt.new, max.time = max.time, dat = dat,
                              nbins = nbins, ndata.types = ndata.types)

  #get marginal loglikel
  pold<- log_marg_likel(alpha = alpha, summary.stats = stats.old,
                      nbins = nbins, ndata.types = ndata.types)
  pnew<- log_marg_likel(alpha = alpha, summary.stats = stats.new,
                      nbins = nbins, ndata.types = ndata.types) + log(p0)
  prob<- exp(pnew - pold)
  rand2<- runif(1)

  if (rand2<prob) return(list(breakpt.new, pnew))
  return(list(breakpt.old, pold))
}
