
#' Summarize observations within bins per track segment
#'
#' Prepares the data that has already been segmented for clustering by Latent
#' Dirichlet Allocation. This function summarizes the counts observed per
#' movement variable bin within each track segment per animal ID.
#'
#' @param dat A data frame of \strong{only} the animal ID, track segment number,
#'   and the discretized data for each movement variable. Animal ID and time
#'   segment must be the first two columns of this data frame. This should be a
#'   simplified form of the output from \code{\link{assign_tseg}}.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{dat}.
#'
#' @return A new data frame that contains the animal ID, track segment number,
#'   and the counts per bin for each movement variable. The names for each of
#'   these bins are labeled according to the order in which the variables were
#'   provided to \code{summarize_tsegs}.
#'
#'
#' @examples
#' #load data
#' data(tracks.seg)
#'
#' #select only id, tseg, SL, and TA columns
#' tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]
#'
#' #run function
#' obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))
#'
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
          dplyr::filter(tseg == k) %>%
          dplyr::select(j+2) %>%
          table()
        mat.list[[j]][k,as.numeric(names(tmp))]<- tmp
      }
    }

    id<- rep(unique(dat.ind$id) %>%
               as.character(), ntseg)
    tseg<- 1:ntseg
    behav.res<- do.call(cbind.data.frame, mat.list) %>%
      data.frame()
    behav.res<- cbind(id, tseg, behav.res)
    behav.res$id<- as.character(behav.res$id)
    obs.list[[i]]<- behav.res
  }
  obs<- dplyr::bind_rows(obs.list)
  obs[is.na(obs)]<- 0  #replace NAs w/ zero
  obs
}
#------------------------------------------------

#' Extract behavior proportion estimates for each track segment
#'
#' Calculates the mean of the posterior for the proportions of each behavior
#' within track segments. These results can be explored to determine the optimal
#' number of latent behavioral states.
#'
#' @param res A list of results returned by \code{\link{cluster_segments}}.
#'   Element \code{theta} stores estimate for behavior proportions for all time
#'   segments.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
#' @param nburn numeric. The length of the burn-in phase.
#' @param nmaxclust numeric. A single number indicating the maximum number of
#'   clusters to test.
#'
#' @return A matrix that stores the proportions of each state/cluster (columns)
#'   per track segment (rows).
#'
#'
#' @examples
#'
#' \donttest{
#' #load data
#' data(tracks.seg)
#'
#' #select only id, tseg, SL, and TA columns
#' tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]
#'
#' #summarize data by track segment
#' obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))
#'
#' #cluster data with LDA
#' res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 1000,
#'                        nburn = 500, nmaxclust = 7, ndata.types = 2)
#'
#' #Extract proportions of behaviors per track segment
#' theta.estim<- extract_prop(res = res, ngibbs = 1000, nburn = 500, nmaxclust = 7)
#' }
#'
#'
#' @export
extract_prop=function(res, ngibbs, nburn, nmaxclust) {
  #Extract and plot proportions of behaviors per track segment
  theta.post<- res$theta[(nburn+1):ngibbs,]  #extract samples from posterior
  theta.estim<- colMeans(theta.post)
  theta.estim<- matrix(data = theta.estim, ncol = nmaxclust) #calc mean of posterior

  theta.estim
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
#' \donttest{
#' #load data
#' data(tracks.seg)
#'
#' #select only id, tseg, SL, and TA columns
#' tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]
#'
#' #summarize data by track segment
#' obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))
#'
#' #cluster data with LDA
#' res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 1000,
#'                        nburn = 500, nmaxclust = 7, ndata.types = 2)
#'
#' #Extract proportions of behaviors per track segment
#' theta.estim<- extract_prop(res = res, ngibbs = 1000, nburn = 500, nmaxclust = 7)
#'
#' #run function
#' behav.res<- get_behav_hist(dat = res, nburn = 500, ngibbs = 1000, nmaxclust = 7,
#'                            var.names = c("Step Length","Turning Angle"))
#' }
#'
#' @importFrom rlang .data
#' @export
get_behav_hist=function(dat, nburn, ngibbs, nmaxclust, var.names) {

  #summarize cluster results by frequency and proportion
  behav.list<- list()
  for (i in 1:length(dat$phi)) {
    tmp<- matrix(dat$phi[[i]][(nburn+1):ngibbs,], length((nburn+1):ngibbs),
                 ncol(dat$phi[[i]]))
    tmp1<- matrix(colMeans(tmp), ncol(tmp) / nmaxclust, nmaxclust, byrow = T)

    behav.list[[i]]<- data.frame(bin = 1:nrow(tmp1), tmp1) %>%
      dplyr::rename_at(dplyr::vars(tidyr::starts_with('X')), ~as.character(1:ncol(tmp1))) %>%
      tidyr::pivot_longer(-.data$bin, names_to = "behav", values_to = "prop") %>%
      dplyr::arrange(.data$behav) %>%
      dplyr::mutate(var = var.names[i])
  }

  #combine params
  behav.res<- dplyr::bind_rows(behav.list)

  behav.res
}
#------------------------------------------------

#' Expand behavior estimates from track segments to observations
#'
#' @param dat A data frame of the animal ID, track segment labels, and all other
#'   data per observation. Animal ID, date, track segment, and observation
#'   number columns must be labeled \emph{id}, \emph{date}, \emph{tseg}, and
#'   \emph{time1}, respectively.
#' @param theta.estim A matrix (returned by \code{\link{extract_prop}})
#'   containing the proportions of each behavioral state as separate columns for
#'   each track segment (rows).
#' @param obs A data frame summarizing the number of observations within each
#'   bin per movement variable that is returned by
#'   \code{\link{summarize_tsegs}}.
#' @param nbehav numeric. The number of behavioral states that will be retained
#'   in 1 to nmaxclust.
#' @param behav.names character. A vector of names to label each state (in
#'   order).
#' @param behav.order numeric. A vector that identifies the order in which the
#'   user would like to rearrange the behavioral states. If satisfied with order
#'   returned by the LDA model, this still must be specified.
#'
#' @return A new data frame that expands behavior proportions for each
#'   observation within all track segments, including the columns labeled
#'   \emph{time1} and \emph{date} from the original \code{dat} data frame.
#'
#'
#' @examples
#'
#' \donttest{
#' #load data
#' data(tracks.seg)
#'
#' #select only id, tseg, SL, and TA columns
#' tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]
#'
#' #summarize data by track segment
#' obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))
#'
#' #cluster data with LDA
#' res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 1000,
#'                        nburn = 500, nmaxclust = 7, ndata.types = 2)
#'
#' #Extract proportions of behaviors per track segment
#' theta.estim<- extract_prop(res = res, ngibbs = 1000, nburn = 500, nmaxclust = 7)
#'
#' #Create augmented matrix by replicating rows (tsegs) according to obs per tseg
#' theta.estim.long<- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
#'                                nbehav = 3, behav.names = c("Encamped","ARS","Transit"),
#'                                behav.order = c(1,2,3))
#' }
#'
#'
#' @importFrom rlang .data
#' @export
expand_behavior=function(dat, theta.estim, obs, nbehav, behav.names, behav.order) {

  #Assign behaviors (via theta) to each track segment
  theta.estim1<- apply(theta.estim[,1:nbehav], 1, function(x) x/sum(x)) %>%
    t()  #normalize probs for only first 3 behaviors being used
  theta.estim1<- data.frame(id = obs$id, tseg = obs$tseg, theta.estim1)
  theta.estim1$id<- as.character(theta.estim1$id)
  names(theta.estim1)<- c("id", "tseg", behav.names)  #define behaviors
  nobs<- data.frame(id = obs$id, tseg = obs$tseg,
                    n = dat %>%
                      dplyr::group_by(.data$id, .data$tseg) %>%
                      dplyr::tally() %>%
                      dplyr::ungroup() %>%
                      dplyr::pull(.data$n))


  for (i in 1:nrow(theta.estim1)) {
    ind<- which(dat$id == theta.estim1$id[i] & dat$tseg == theta.estim1$tseg[i])

    if (i == 1) {
      theta.estim2<- rep(theta.estim1[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim1), byrow = TRUE)
    } else {
      tmp<- rep(theta.estim1[i,], nobs$n[i]) %>%
        unlist() %>%
        matrix(nrow = nobs$n[i], ncol = ncol(theta.estim1), byrow = TRUE)

      theta.estim2<- rbind(theta.estim2, tmp)
    }
  }

  colnames(theta.estim2)<- names(theta.estim1)
  theta.estim2<- data.frame(theta.estim2, time1 = dat$time1, date = dat$date,
                            stringsAsFactors = FALSE)

  #Change col classes to numeric besides ID
  ind1<- which(names(theta.estim2) != "id")
  theta.estim2<- theta.estim2 %>%
    dplyr::mutate_at(names(theta.estim2)[ind1], as.numeric) %>%
    dplyr::select(.data$id, .data$tseg, .data$time1, .data$date, dplyr::everything())

  #Change into long format
  theta.estim.long<- tidyr::pivot_longer(theta.estim2, cols = -c(1:4),
                                         names_to = "behavior", values_to = "prop")
  theta.estim.long$behavior<- factor(theta.estim.long$behavior,
                                     levels = behav.names[behav.order])
  theta.estim.long<- theta.estim.long %>%
    dplyr::arrange(.data$behavior) %>%
    dplyr::mutate_at("date", lubridate::as_datetime)

  theta.estim.long
}
#------------------------------------------------

#' Assign behavior estimates to observations
#'
#' @param dat.orig A data frame that contains all of the original data for all
#'   animal IDs. Must be same as was used to originally segment the tracks. Must
#'   have columns \code{obs} and \code{time1} generated by
#'   \code{\link{filter_time}}.
#' @param dat.seg.list A list of data associated with each animal ID where names
#'   of list elements are the ID names and tracks have already been segmented.
#'   Must have columns \code{obs} and \code{time1} generated by
#'   \code{\link{filter_time}}.
#' @param theta.estim.long A data frame in long format where each observation
#'   (\emph{time1}) of each track segment (\emph{tseg}) of each animal ID
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
#'   given track segment for which the observation belongs (\code{behav}). This
#'   is merged with the original data frame \code{dat.orig}, so any observations
#'   that were excluded (not at primary time interval) will show \code{NA} for
#'   behavior estimates.
#'
#'
#' @examples
#'
#' \donttest{
#' #load original and segmented data
#' data(tracks)
#' data(tracks.seg)
#'
#' #convert segmented dataset into list
#' tracks.list<- df_to_list(dat = tracks.seg, ind = "id")
#'
#' #select only id, tseg, SL, and TA columns
#' tracks.seg2<- tracks.seg[,c("id","tseg","SL","TA")]
#'
#' #summarize data by track segment
#' obs<- summarize_tsegs(dat = tracks.seg2, nbins = c(5,8))
#'
#' #cluster data with LDA
#' res<- cluster_segments(dat = obs, gamma1 = 0.1, alpha = 0.1, ngibbs = 1000,
#'                        nburn = 500, nmaxclust = 7, ndata.types = 2)
#'
#' #Extract proportions of behaviors per track segment
#' theta.estim<- extract_prop(res = res, ngibbs = 1000, nburn = 500, nmaxclust = 7)
#'
#' #Create augmented matrix by replicating rows (tsegs) according to obs per tseg
#' theta.estim.long<- expand_behavior(dat = tracks.seg, theta.estim = theta.estim, obs = obs,
#'                                nbehav = 3, behav.names = c("Encamped","ARS","Transit"),
#'                                behav.order = c(1,2,3))
#'
#' #Run function
#' dat.out<- assign_behavior(dat.orig = tracks, dat.seg.list = tracks.list,
#'                           theta.estim.long = theta.estim.long,
#'                           behav.names = c("Encamped","ARS","Transit"))
#' }
#'
#'
#' @export
assign_behavior=function(dat.orig, dat.seg.list, theta.estim.long, behav.names) {

  for (i in 1:length(dat.seg.list)) {
    sub<- theta.estim.long[theta.estim.long$id == unique(dat.seg.list[[i]]$id),]
    sub<- sub %>%
      dplyr::arrange(.data$tseg, .data$date, .data$behavior) %>%
      tidyr::pivot_wider(names_from = .data$behavior, values_from = .data$prop)
    sub<- sub %>%
      dplyr::mutate(behav = behav.names[apply(sub[,5:ncol(sub)], 1, which.max)])


    dat.seg.list[[i]]<- cbind(dat.seg.list[[i]], sub[,5:ncol(sub)])
  }

  #Convert original data to list and filter for only IDs analyzed
  dat.orig.list<- df_to_list(dat = dat.orig, ind = "id")
  ind<- which(names(dat.orig.list) %in% names(dat.seg.list))
  dat.orig.list<- dat.orig.list[ind]

  #Add obs col to dat.orig.list
  dat.orig.list1<- purrr::map(dat.orig.list, ~dplyr::mutate(.x, obs = 1:nrow(.x)))

  #Merge behav estimates w/ original data
  dat1<- purrr::map2(dat.orig.list1, dat.seg.list,
                     ~dplyr::left_join(.x, .y[,c("obs","tseg", behav.names, "behav")],
                                       by = "obs")
  ) %>%
    dplyr::bind_rows()

  dat1$behav<- factor(dat1$behav, levels = behav.names)
  dat1
}
