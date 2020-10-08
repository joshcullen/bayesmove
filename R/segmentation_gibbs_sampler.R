#' Internal function that runs RJMCMC on a single animal ID
#'
#' This function serves as a wrapper for \code{\link{samp_move}} by running this
#' sampler for each iteration of the MCMC chain. It is called by
#' \code{\link{segment_behavior}} to run the RJMCMC on all animal IDs
#' simultaneously.
#'
#' @param dat A data frame that only contains columns for the animal IDs and for
#'   each of the discretized movement variables.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{dat}.
#' @param alpha numeric. A single value used to specify the hyperparameter for
#'   the prior distribution. A standard value for \code{alpha} is typically 1,
#'   which corresponds with a vague prior on the Dirichlet distribution.
#' @param breakpt numeric. A vector of breakpoints if pre-specifying where they
#'   may occur, otherwise \code{NULL}.
#'
#' @return A list of the breakpoints, the number of breakpoints, and the log
#'   marginal likelihood at each MCMC iteration, as well as the time it took the
#'   model to finish running. This is only provided for the data of a single
#'   animal ID.
#'
#'
#'
#'
#'
behav_gibbs_sampler=function(dat, ngibbs, nbins, alpha, breakpt) {

  start.time<- Sys.time()  #start timer

  uni.id<- unique(dat$id)  #need this query if running multiple IDs in parallel
  dat<- subset(dat, select = -id)

  #useful stuff
  max.time<- nrow(dat)
  ndata.types<- length(nbins)

  #starting values
  if (is.null(breakpt)) breakpt<- floor(max.time/2)

  #to store results
  res.brks<- vector("list", ngibbs)
  res.LML<- matrix(NA, 1, (ngibbs+1))
  res.nbrks<- matrix(NA, 1, (ngibbs+1))
  store.param<- matrix(NA, ngibbs, 2)

  for (i in 1:ngibbs){
    vals<- samp_move(breakpt = breakpt, max.time = max.time, dat = dat,
                   alpha = alpha, nbins = nbins, ndata.types = ndata.types)

    breakpt<- vals[[1]]

    #store results
    res.brks[[i]]<- breakpt
    store.param[i,]<- c(length(breakpt), vals[[2]])

  }

  tmp<- store.param[,1]
  res.nbrks[1,]<- c(uni.id,tmp)
  colnames(res.nbrks)<- c('id', paste0("Iter_",1:ngibbs))

  tmp<- store.param[,2]
  res.LML[1,]<- c(uni.id,tmp)
  colnames(res.LML)<- c('id', paste0("Iter_",1:ngibbs))

  end.time<- Sys.time()
  elapsed.time<- difftime(end.time, start.time, units = "min")  #end timer

  list(breakpt=res.brks, nbrks=res.nbrks, LML=res.LML, elapsed.time=elapsed.time)
}
#----------------------------------------------------

#' Segmentation model to estimate breakpoints
#'
#' This function performs the reversible-jump MCMC algorithm using a Gibbs
#' sampler, which estimates the breakpoints of the movement variables for each
#' of the animal IDs. This is the first stage of the two-stage Bayesian model that
#' estimates proportions of behavioral states by first segmenting individual
#' tracks into relatively homogeneous segments of movement.
#'
#' This model is run in parallel using the \code{future} package. To ensure that
#' the model is run in parallel, the \code{\link[future]{plan}} must be used
#' with \code{future::multisession} as the argument for most operating systems.
#' Otherwise, model will run sequentially by default if this is not set before
#' running \code{segment_behavior}.
#'
#' @param data A list where each element stores the data for a separate animal
#'   ID. List elements are data frames that only contain columns for the animal
#'   ID and for each of the discretized movement variables.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{data}.
#' @param alpha numeric. A single value used to specify the hyperparameter for
#'   the prior distribution. A standard value for \code{alpha} is typically 1,
#'   which corresponds with a vague prior on the Dirichlet distribution.
#' @param breakpt A list where each element stores a vector of breakpoints if
#'   pre-specifying where they may occur for each animal ID. By default this is
#'   set to \code{NULL}.
#'
#' @return A list of model results is returned where elements include the
#'   breakpoints, number of breakpoints, and log marginal likelihood at each
#'   iteration of the MCMC chain for all animal IDs. The time it took the model
#'   to finish running for each animal ID are also stored and returned.
#'
#'
#' @examples
#' \donttest{
#' #load data
#' data(tracks.list)
#'
#' #subset only first track
#' tracks.list<- tracks.list[1]
#'
#' #only retain id and discretized step length (SL) and turning angle (TA) columns
#' tracks.list2<- purrr::map(tracks.list,
#'                    subset,
#'                   select = c(id, SL, TA))
#'
#'
#' set.seed(1)
#'
#' # Define model params
#' alpha<- 1
#' ngibbs<- 1000
#' nbins<- c(5,8)
#'
#' future::plan(future::multisession)  #run all MCMC chains in parallel
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                            alpha = alpha)
#'
#' future::plan(future::sequential)  #return to single core
#' }
#'
#' @importFrom future "plan"
#'
#' @export
segment_behavior=function(data, ngibbs, nbins, alpha,
                          breakpt = purrr::map(names(data), ~ NULL)) {

  tictoc::tic()  #start timer
  mod<- furrr::future_map2(data, breakpt,
                           ~behav_gibbs_sampler(dat = .x, ngibbs = ngibbs, nbins = nbins,
                                                         alpha = alpha, breakpt = .y),
                   .progress = TRUE, .options = furrr::future_options(seed = T))
  tictoc::toc()  #provide elapsed time




  brkpts<- purrr::map(mod, 1)  #create list of all sets breakpoints by ID


  nbrks<- purrr::map_dfr(mod, 2) %>%
    unlist() %>%
    matrix(.data, nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.nbrks<- ncol(nbrks)
  nbrks[,2:ncol.nbrks]<- apply(nbrks[,2:ncol.nbrks], 2, function(x) as.numeric(as.character(x)))


  LML<- purrr::map_dfr(mod, 3) %>%
    unlist() %>%
    matrix(.data, nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.LML<- ncol(LML)
  LML[,2:ncol.LML]<- apply(LML[,2:ncol.LML], 2, function(x) as.numeric(as.character(x)))


  elapsed.time<- purrr::map_dfr(mod, 4) %>%
    t() %>%
    data.frame()  #create DF of elapsed time
  names(elapsed.time)<- "time"
  elapsed.time<- elapsed.time %>%
    dplyr::mutate_at("time", as.character)


  list(brkpts = brkpts, nbrks = nbrks, LML = LML, elapsed.time = elapsed.time)
}
