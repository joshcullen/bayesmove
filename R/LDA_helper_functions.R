
#' Summarize observations within bins per time segment
#'
#' Pepares the data that has already been segmented for clustering by Latent
#' Dirichlet Allocation. This function summarizes the counts observed per
#' movement variable bin within each time segment per animal ID.
#'
#' @param dat A data frame of the animal ID, time segment number, and the
#'   discretized data for each movement variable. Animal ID and time segment
#'   must be the first two columns of this data frame.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{dat}.
#'
#' @return A new data frame that contains the animal ID, time segment number,
#'   and the counts per bin for each movement variable. The names for each of
#'   these bins are labeled according to the order in which the variables were
#'   provided to \code{summarize_tsegs}.
#'
#'
#' @examples
#'
#' @export
summarize_tsegs=function(dat, nbins){

  #create list of input and to store output
  dat.list<- df_to_list(dat = dat, ind = "id")
  id<- unique(dat$id) %>%
    as.character()
  n<- length(id)
  obs.list<- vector("list", n)
  names(obs.list)<- id

  #calculate # of obs in each bin (per move var) by tseg
  for (i in 1:length(dat.list)) {
    dat.ind=dat.list[[i]]
    ntseg=max(dat.ind$tseg)


    #generate counts of obs per var bin by tseg for all vars; 'id' and 'tseg'
    #need to be 1st two cols
    mat.list=list()
    for (j in 1:length(nbins)) {
      mat.list[[j]]<- matrix(0, ntseg, nbins[j])
      colnames(mat.list[[j]])<- paste0("y",j,".",1:nbins[j])
      for (k in 1:ntseg){
        tmp<- dat.ind %>%
          filter(tseg == k) %>%
          dplyr::select(j+2) %>%
          table()
        mat.list[[j]][k,as.numeric(names(tmp))]<- tmp
      }
    }

    id<- rep(unique(dat.ind$id) %>%
               as.character(), ntseg)
    tseg<- 1:ntseg
    behav.res<- do.call(cbind.data.frame, mat.list) %>%
      data.frame() %>%
      cbind(id, tseg, .)
    behav.res$id<- as.character(behav.res$id)
    obs.list[[i]]<- behav.res
  }
  obs<- bind_rows(obs.list)
  obs[is.na(obs)]<- 0  #replace NAs w/ zero
  obs
}

#------------------------------------------------

#' Extract bin estimates from Latent Dirichlet Allocation model
#'
#' Pulls model results for the estimates of bin proportions per movement
#' variable from the posterior distribution. This can be used for visualization
#' of movement variable distribution for each behavior estimated.
#'
#' @param dat The list object returned by the LDA model
#'   (\code{\link{cluster_segments}}). Used for extracting the element
#'   \emph{phi}.
#' @param nburn numeric. The length of the burn-in phase.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
#' @param nmaxclust numeric. The maximum number of clusters on which to
#'   attribute behaviors.
#' @param var.names character. A vector of names used for each of the movement
#'   variables. Must be in the same order as were listed within the data frame
#'   returned by \code{\link{summarize_tsegs}}.
#'
#' @return A data frame that contains columns for bin number, behavioral state,
#'   proportion represented by a given bin, and movement variable name. This is
#'   displayed in a long format, which is easier to visualize using
#'   \code{ggplot2}.
#'
#'
#' @examples
#'
#' @export
get_behav_hist=function(dat, nburn, ngibbs, nmaxclust, var.names) {

  #summarize cluster results by frequency and proportion
  behav.list<- list()
  for (i in 1:length(dat$phi)) {
    tmp<- matrix(dat$phi[[i]][(nburn+1):ngibbs,], length((nburn+1):ngibbs), ncol(dat$phi[[i]]))
    tmp1<- matrix(colMeans(tmp), ncol(tmp) / nmaxclust, nmaxclust, byrow = T)
    behav.list[[i]]<- data.frame(bin = 1:nrow(tmp1), tmp1) %>%
      rename_at(vars(starts_with('X')), ~as.character(1:ncol(tmp1))) %>%
      pivot_longer(-bin, names_to = "behav", values_to = "prop") %>%
      arrange(behav) %>%
      mutate(var = var.names[i])
  }

  #combine params
  behav.res<- bind_rows(behav.list)

  behav.res
}
#------------------------------------------------

#' Expand behavior estimates from time segments to observations
#'
#' @param dat A data frame of the animal ID, time segment labels, and all other
#'   data per observation. Animal ID and time segment columns must be labeled
#'   \emph{id} and \emph{tseg}, respectively.
#' @param theta.estim A data frame containing the animal ID, time segment, and
#'   proportions of each behavioral state as separate columns. Animal ID and
#'   time segment columns must be labeled \emph{id} and \emph{tseg},
#'   respectively.
#' @param nobs A data frame containing the animal ID, time segment, and number
#'   of observations recorded within a given time segment. Animal ID, time
#'   segment, and sample size columns must be labeled \emph{id}, \emph{tseg},
#'   and \emph{n}, respectively.
#'
#' @return A new data frame that expands behavior proportions for each
#'   observation within all time segments, including the columns labeled
#'   \emph{time1} and \emph{date} from the original \code{dat} data frame.
#'
#'
#' @examples
#'
#' @export
expand_behavior=function(dat, theta.estim, nobs) {

  for (i in 1:nrow(theta.estim)) {
    ind<- which(dat$id == theta.estim$id[i] & dat$tseg == theta.estim$tseg[i])

    if (i == 1) {
      theta.estim2<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)
    } else {
      tmp<- rep(theta.estim[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim), byrow = TRUE)

      theta.estim2<- rbind(theta.estim2, tmp)
    }
  }

  colnames(theta.estim2)<- names(theta.estim)
  theta.estim2<- data.frame(theta.estim2, time1 = dat$time1, date = dat$date,
                            stringsAsFactors = FALSE)

  theta.estim2
}
#------------------------------------------------

#' Assign behavior estimates to observations
#'
#' @param dat.list A list of data associated with each animal ID where names of
#'   list elements are the ID names.
#' @param theta.estim.long A data frame in long format where each observation
#'   (\emph{time1}) of each time segment (\emph{tseg}) of each animal ID
#'   (\emph{id}) has separate rows for behavior proportion estimates per state.
#'   Columns for behavior and proportion estimates should be labeled
#'   \emph{behavior} and \emph{prop}, respectively. Date (in POSIXct format)
#'   should also be included as a column labeled \emph{date}.
#' @param behav.names character. A vector of names to label each state (in
#'   order).
#'
#' @return A data frame of all animal IDs where columns (with names from
#'   \code{behav.names}) include proportions of each behavioral state per
#'   observation, as well as a column that stores the dominant behavior within a
#'   given time segment for which the observation belongs (\code{behav}).
#'
#'
#' @examples
#'
#' @export
assign_behavior=function(dat.list, theta.estim.long, behav.names) {  #assign dominant behavior to observations

  for (i in 1:length(dat.list)) {
    sub<- theta.estim.long[theta.estim.long$id == unique(dat.list[[i]]$id),]
    sub<- sub %>%
      arrange(tseg, date, behavior) %>%
      pivot_wider(names_from = behavior, values_from = prop) %>%
      mutate(behav = behav.names[apply(.[,5:ncol(.)], 1, which.max)])


    dat.list[[i]]<- cbind(dat.list[[i]], sub[,5:ncol(sub)])
  }

  #Convert to DF
  dat2<- dplyr::bind_rows(dat.list)

  dat2
}
