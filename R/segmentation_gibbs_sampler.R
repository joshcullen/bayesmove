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
#' @export
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
#' with \code{multisession} as the argument for most operating systems.
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
#' \dontrun{
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
#' date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
#' dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
#' dt<- c(dt, NA)
#' var<- rep(sample(c(1,2), 100, replace = TRUE), each = 10)
#' id<- rep(1:10, each = 100)
#'
#'
#' #create data frame and round time
#' dat<- data.frame(id, date, dt, step, angle, var)
#' dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC")
#'
#'
#' #define limits for each bin
#' dist.lims<- c(quantile(step, c(0, 0.25, 0.5, 0.75, 0.95)), max(step))  #5 bins
#' angle.lims<- c(-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi)  #8 bins
#'
#' #discretize step and angle
#' dat1<- discrete_move_var(dat = dat, lims = list(dist.lims, angle.lims),
#'                          varIn = c("step", "angle"),
#'                          varOut = c("SL","TA"))
#'
#'
#' #create list and filter by primary time step
#' dat.list<- df_to_list(dat = dat1, ind = "id")
#' dat.list.filt<- filter_time(dat.list = dat.list, int = 3600)
#'
#'
#' #find breakpoints to pre-specify by ID
#' breaks<- lapply(dat.list.filt, find_breaks, ind = "var")
#'
#'
#'
#'
#'
#'
#' #perform segmentation w/o pre-specifying breakpoints
#' dat.list.filt1<- lapply(dat.list.filt,
#'                         function(x) subset(x, select = c(id, SL, TA)))
#' future::plan(future::multisession)
#' dat.res1<- segment_behavior(data = dat.list.filt1, ngibbs = 1000, nbins = c(5,8),
#'                             alpha = 1)
#'
#'
#' #perform segmentation w/ pre-specifying breakpoints
#' dat.list.filt2<- lapply(dat.list.filt,
#'                         function(x) subset(x, select = c(id, SL, TA, var)))
#' future::plan(future::multisession)
#' dat.res2<- segment_behavior(data = dat.list.filt2, ngibbs = 1000, nbins = c(5,8,2),
#'                             alpha = 1, breakpt = breaks)
#'
#' future::plan(future::sequential)
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
                   .progress = TRUE)
  tictoc::toc()  #provide elapsed time




  brkpts<- purrr::map(mod, 1)  #create list of all sets breakpoints by ID


  nbrks<- purrr::map_dfr(mod, 2) %>%
    unlist() %>%
    matrix(., nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of number of breakpoints by ID
  names(nbrks)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.nbrks<- ncol(nbrks)
  nbrks<- nbrks %>%
    dplyr::mutate_at(2:ncol.nbrks, as.character) %>%
    dplyr::mutate_at(2:ncol.nbrks, as.numeric) %>%
    dplyr::mutate_at(1, as.character)


  LML<- purrr::map_dfr(mod, 3) %>%
    unlist() %>%
    matrix(., nrow = length(mod), ncol = (ngibbs + 1), byrow = T) %>%
    data.frame()  #create DF of LML by ID
  names(LML)<- c('id', paste0("Iter_", 1:ngibbs))
  ncol.LML<- ncol(LML)
  LML<- LML %>%
    dplyr::mutate_at(2:ncol.LML, as.character) %>%
    dplyr::mutate_at(2:ncol.LML, as.numeric) %>%
    dplyr::mutate_at(1, as.character)


  elapsed.time<- purrr::map_dfr(mod, 4) %>%
    t() %>%
    data.frame()  #create DF of elapsed time
  names(elapsed.time)<- "time"
  elapsed.time<- elapsed.time %>%
    dplyr::mutate_at("time", as.character)


  list(brkpts = brkpts, nbrks = nbrks, LML = LML, elapsed.time = elapsed.time)
}
