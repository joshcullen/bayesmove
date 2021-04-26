#' Cluster observations into behavioral states
#'
#' This function uses a Gibbs sampler within a mixture model to estimate the
#' optimal number of behavioral states, the state-dependent distributions, and
#' to assign behavioral states to each observation. This model does not assume
#' an underlying mechanistic process.
#'
#' The mixture model analyzes all animal IDs pooled together, thus providing a
#' population-level estimate of behavioral states.
#'
#' @param dat A data frame that **only** contains columns for the discretized
#'   movement variables.
#' @param alpha numeric. A single value used to specify the hyperparameter for
#'   the prior distribution.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#' @param nburn numeric. The length of the burn-in phase.
#'
#' @return A list of model results is returned where elements include the
#'   \code{phi} matrix for each data stream, \code{theta} matrix, log likelihood
#'   estimates for each iteration of the MCMC chain \code{loglikel}, a list of
#'   the MAP estimates of the latent states for each observation \code{z.MAP}, a
#'   matrix of the whole posterior of state assignments per observation
#'   \code{z.posterior}, and a vector \code{gamma1} of estimates for the gamma
#'   hyperparameter.
#' @examples
#' \donttest{
#' data(tracks.list)
#'
#' #convert from list to data frame
#' tracks.list<- dplyr::bind_rows(tracks.list)
#'
#' #only retain id and discretized step length (SL) and turning angle (TA) columns
#' tracks<- subset(tracks.list, select = c(SL, TA))
#'
#'
#' set.seed(1)
#'
#' # Define model params
#' alpha=0.1
#' ngibbs=1000
#' nburn=ngibbs/2
#' nmaxclust=7
#'
#' dat.res<- cluster_obs(dat = tracks, alpha = alpha, ngibbs = ngibbs,
#'                            nmaxclust = nmaxclust, nburn = nburn)
#' }
#'
#' @export
cluster_obs=function(dat, alpha, ngibbs, nmaxclust, nburn){

  nobs=nrow(dat)
  ndata.types=ncol(dat)


  #initial values
  z=sample(1:nmaxclust,size=nobs,replace=T)
  nbins=apply(dat,2,max,na.rm=T)
  phi=list()
  for (i in 1:ndata.types){
    phi[[i]]=matrix(1/nbins[i],nmaxclust,nbins[i])
  }
  theta=rep(1/nmaxclust,nmaxclust)
  gamma1=0.1

  #get nmat
  nmat=list()
  for (i in 1:ndata.types){
    nmat[[i]]=SummarizeDat(z=z-1, dat=dat[,i]-1, ncateg=nbins[i],nbehav=nmaxclust, nobs=nobs)
  }

  #prepare for gibbs
  store.phi=list()
  for (i in 1:ndata.types){
    store.phi[[i]]=matrix(NA,ngibbs,nmaxclust*nbins[i])
  }
  store.theta=matrix(NA,ngibbs,nmaxclust)
  store.loglikel=rep(NA,ngibbs)
  store.gamma1=rep(NA,ngibbs)
  store.z=matrix(0,nobs,nmaxclust)

  #run gibbs sampler
  max.llk=-Inf
  gamma.possib=seq(from=0.1,to=1,by=0.05) #possible values for gamma

  #progress bar
  pb <- progress::progress_bar$new(
    format = " iteration (:current/:total) [:bar] :percent [Elapsed: :elapsed, Remaining: :eta]",
    total = ngibbs, clear = FALSE, width = 100)


  for (i in 1:ngibbs){
    pb$tick()  #create progress bar

    #sample from FCD's
    lphi=list()
    for (j in 1:ndata.types) lphi[[j]]=log(phi[[j]])
    ltheta=log(theta)
    z=sample.z.mixmod(nobs=nobs,nmaxclust=nmaxclust,dat=dat,ltheta=ltheta,lphi=lphi,
               ndata.types=ndata.types)
    for (j in 1:ndata.types){
      nmat[[j]]=SummarizeDat(z=z-1, dat=dat[,j]-1, ncateg=nbins[j],nbehav=nmaxclust,
                             nobs=nobs)
    }

    tmp=sample.v.mixmod(z=z,gamma1=gamma1,nmaxclust=nmaxclust)
    theta=tmp$theta
    v=tmp$v

    phi=sample.phi.mixmod(alpha=alpha,nmaxclust=nmaxclust,
                   nbins=nbins,ndata.types=ndata.types,nmat=nmat)

    gamma1=sample.gamma.mixmod(v=v,ngroup=nmaxclust,gamma.possib=gamma.possib)

    #calculate log-likelihood
    llk=get.llk.mixmod(phi=phi,theta=theta,ndata.types=ndata.types,dat=dat,
                nobs=nobs,nmaxclust=nmaxclust)

    #store results
    for (j in 1:ndata.types){
      store.phi[[j]][i,]=phi[[j]]
    }
    store.theta[i,]=theta
    store.loglikel[i]=llk
    store.gamma1[i]=gamma1

    if (i> nburn){
      store.z=StoreZ(z=z-1,store_z=store.z,nobs=nobs)
    }

    #re-order clusters
    if (i < nburn & i%%50==0){
      ordem=order(theta,decreasing=T)
      theta=theta[ordem]

      for (j in 1:ndata.types){
        phi[[j]]=phi[[j]][ordem,]
        nmat[[j]]=nmat[[j]][ordem,]
      }

      znew=z
      for (j in 1:nmaxclust){
        cond=z==ordem[j]
        znew[cond]=j
      }
      z=znew
    }
  }

  if (i > nburn & llk>max.llk){
    z.MAP=z
    max.llk=llk
  }

  list(phi=store.phi,theta=store.theta,
       loglikel=store.loglikel,z.MAP=z.MAP,
       z.posterior=store.z,
       gamma1=store.gamma1)
}
