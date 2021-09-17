#' Discretize movement variables
#'
#' Convert movement variables from continuous to discrete values for analysis by
#' \code{\link{segment_behavior}}.
#'
#' @param dat A data frame that contains the variable(s) of interest to convert
#'   from continuous to discrete values.
#' @param lims A list of the bin limits for each variable. Each element of the
#'   list should be a vector of real numbers.
#' @param varIn A vector of names for the continuous variable stored as columns
#'   within \code{dat}.
#' @param varOut A vector of names for the storage of the discrete variables
#'   returned by the function.
#'
#' @return A data frame with new columns of discretized variables as labeled by
#'   \code{varOut}.
#'
#' @import dplyr
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #subset only first track
#' tracks<- tracks[tracks$id == "id1",]
#'
#' #calculate step lengths and turning angles
#' tracks<- prep_data(dat = tracks, coord.names = c("x","y"), id = "id")
#'
#' #round times to nearest interval of interest (e.g. 3600 s or 1 hr)
#' tracks<- round_track_time(dat = tracks, id = "id", int = 3600, tol = 180, time.zone = "UTC",
#'                                  units = "secs")
#'
#' #create list from data frame
#' tracks.list<- df_to_list(dat = tracks, ind = "id")
#'
#' #filter observations to only 1 hr (or 3600 s)
#' tracks_filt.list<- filter_time(dat.list = tracks.list, int = 3600)
#'
#' #define bin number and limits for turning angles and step lengths
#' angle.bin.lims=seq(from=-pi, to=pi, by=pi/4)  #8 bins
#' dist.bin.lims=quantile(tracks[tracks$dt == 3600,]$step,
#'                       c(0,0.25,0.50,0.75,0.90,1), na.rm=TRUE)  #5 bins
#'
#'
#' # Assign bins to observations
#' tracks_disc.list<- purrr::map(tracks_filt.list,
#'                       discrete_move_var,
#'                       lims = list(dist.bin.lims, angle.bin.lims),
#'                       varIn = c("step", "angle"),
#'                       varOut = c("SL", "TA"))
#'
#'
#' @export
discrete_move_var=function(dat, lims, varIn, varOut){

  for (i in 1:length(lims)) {
    for(j in 1:length(lims[[i]])) {
      tmp<- which(dat[,varIn[i]] >= lims[[i]][j] & dat[,varIn[i]] < lims[[i]][j+1])
      dat[tmp,varOut[i]]<- j
    }
    tmp<- which(dat[,varIn[i]] == lims[[i]][length(lims[[i]])])
    dat[tmp,varOut[i]]<- length(lims[[i]]) - 1
  }

  dat
}

#---------------------------------------------


#' Round time to nearest interval
#'
#' Rounds sampling intervals that are close, but not exactly the time interval
#' of interest (e.g., 240 s instead of 300 s). This can be performed on multiple
#' time intervals, but only using a single tolerance value. This function
#' prepares the data to be analyzed by \code{\link{segment_behavior}}, which
#' requires that all time intervals exactly match the primary time interval when
#' analyzing step lengths and turning angles. Columns storing the time intervals
#' and dates must be labeled \code{dt} and \code{date}, respectively, where
#' dates are of class \code{POSIXct}.
#'
#' @param dat A data frame that contains the sampling interval of the
#'   observations.
#' @param id character. The name of the column storing the animal IDs.
#' @param int numeric. A vector of the time interval(s) of on which to perform
#'   rounding.
#' @param tol numeric. A single tolerance value on which to round any \code{int}
#'   that were specified.
#' @param time.zone character. Specify the time zone for which the date-times
#'   were recorded. Set to UTC by default. Refer to \code{base::OlsonNames} to view
#'   all possible time zones.
#' @param units character. The units of the selected time interval \code{int},
#'   which can be selected from one of "secs", "mins", "hours", "days", or
#'   "weeks".
#'
#' @return A data frame where \code{dt} and \code{date} are both adjusted based
#'   upon the rounding of time intervals according to the specified tolerance.
#'
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #subset only first track
#' tracks<- tracks[tracks$id == "id1",]
#'
#' #calculate step lengths and turning angles
#' tracks<- prep_data(dat = tracks, coord.names = c("x","y"), id = "id")
#'
#' #round times to nearest interval of interest (e.g. 3600 s or 1 hr)
#' tracks<- round_track_time(dat = tracks, id = "id", int = 3600, tol = 180, time.zone = "UTC",
#'                           units = "secs")
#'
#' @export
round_track_time = function(dat, id, int, tol, time.zone = "UTC", units) {

  dat<- df_to_list(dat, ind = id)
  for (i in 1:length(dat)) {
    tmp<- matrix(NA, nrow(dat[[i]]), 1)

    if (length(int) == 1) {  #when using only 1 time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j,]<- NA
        } else if (dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                   dat[[i]]$dt[j] != int) {
          tmp[j,]<- int
        } else {
          tmp[j,]<- dat[[i]]$dt[j]
        }
      }
    } else {  #when using more than one time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j,]<- NA
        } else if (sum(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                       dat[[i]]$dt[j] != int) > 0) {
          ind<- which(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                        dat[[i]]$dt[j] != int)
          tmp[j,]<- int[ind]
        } else {
          tmp[j,]<- dat[[i]]$dt[j]
        }
      }
    }
    dat[[i]]$dt<- tmp[,1]


    if (units == "secs") {
      tmp.dt<- lubridate::seconds(dat[[i]]$dt)
    } else if (units == "mins") {
      tmp.dt<- lubridate::seconds(dat[[i]]$dt)*60
    } else if (units == "hours") {
      tmp.dt<- lubridate::seconds(dat[[i]]$dt)*60*60
    } else if (units == "days") {
      tmp.dt<- lubridate::seconds(dat[[i]]$dt)*60*60*24
    } else if (units == "weeks") {
      tmp.dt<- lubridate::seconds(dat[[i]]$dt)*60*60*24*7
    } else {
      stop("Units must be either secs, mins, hours, days, or weeks")
    }
    tmp.date<- cumsum(c(as.numeric(dat[[i]]$date[1]), tmp.dt[-length(tmp.dt)]))
    dat[[i]]$date<- tmp.date %>%
      as.POSIXct(origin = '1970-01-01', tz = time.zone)
  }

  dat<- dplyr::bind_rows(dat)
  dat
}
#---------------------------------------------

#' Filter observations for time interval of interest
#'
#' Selects observations that belong to the time interval of interest and removes
#' all others. This function also removes entire IDs from the dataset when there
#' is one or fewer observations at this time interval. This function works
#' closely with \code{\link{round_track_time}} to only retain observations
#' sampled at a regular time interval, which is important for analyzing step
#' lengths and turning angles. Column storing the time intervals must be labeled
#' \code{dt}.
#'
#' @param dat.list A list of data associated with each animal ID where names of
#'   list elements are the ID names.
#' @param int numeric. The time interval of interest.
#'
#' @return A list where observations for each animal ID (element) has been
#'   filtered for \code{int}. Two columns (\code{obs} and \code{time1}) are
#'   added for each list element (ID), which store the original observation
#'   number before filtering and the new observation number after filtering,
#'   respectively.
#'
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #subset only first track
#' tracks<- tracks[tracks$id == "id1",]
#'
#' #calculate step lengths and turning angles
#' tracks<- prep_data(dat = tracks, coord.names = c("x","y"), id = "id")
#'
#' #round times to nearest interval of interest (e.g. 3600 s or 1 hr)
#' tracks<- round_track_time(dat = tracks, id = "id", int = 3600, tol = 180, time.zone = "UTC",
#'                               units = "secs")
#'
#' #create list from data frame
#' tracks.list<- df_to_list(dat = tracks, ind = "id")
#'
#' #filter observations to only 1 hr (or 3600 s)
#' tracks_filt.list<- filter_time(dat.list = tracks.list, int = 3600)
#'
#' @importFrom rlang .data
#' @export
filter_time=function(dat.list, int) {

  id<- names(dat.list)

  cond<- matrix(0, length(dat.list), 1)
  for (i in 1:length(dat.list)) {  #don't include IDs w <= 1 obs of dt == int
    cond[i,]<- if(nrow(dat.list[[i]][dat.list[[i]]$dt == int,]) > 1) {
      1
    } else {
      0
    }
  }
  tmp<- cbind(id, cond)
  ind<- which(tmp[,2] == 1)

  n<- length(ind)
  behav.list<- vector("list", n)
  names(behav.list)<- id[ind]

  for (i in 1:n) {
    behav.list[[i]]<- dat.list[[ind[i]]] %>%
      dplyr::mutate(obs = 1:length(.data$id)) %>%
      dplyr::filter(.data$dt == int) %>%
      dplyr::mutate(time1 = 1:length(.data$id))
  }

  behav.list
}
#---------------------------------------------

#' Internal function that transforms a vector of bin numbers to a
#' presence-absence matrix
#'
#' Transforms vectors of bin numbers into full matrices for plotting as a
#' heatmap.
#'
#' @param dat A data frame for a single animal ID that contains only columns for
#'   the ID and each of the movement variables that were analyzed by
#'   \code{\link{segment_behavior}}. The ID column must be first.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{data}.
#'
#' @return A list where each element stores the presence-absence matrix for each
#'   of the movement variables.
#'
#'
#'
#'
behav_seg_image=function(dat, nbins) {

  dat1<- data.frame(dat[,-1])  #remove id col
  names(dat1)<- names(dat)[-1]
  behav.list<- purrr::map2(list(dat1), nbins, ~matrix(0, nrow = nrow(.x), ncol = .y))
  for (i in 1:length(behav.list)) {
    for (j in 1:nrow(dat1)){
      behav.list[[i]][j,dat1[,i][j]]=1
    }
  }

  names(behav.list)<- names(dat1)
  behav.list
}
#---------------------------------------

#' Internal function that adds segment numbers to observations
#'
#' After breakpoints have been extracted for each animal ID, this function
#' assigns the associated segment number to observations for each animal ID.
#' These segments of observations will be used in the second stage of the model
#' framework to perform mixed-membership clustering by Latent Dirichlet
#' Allocation.
#'
#' @param dat A data frame that contains all data associated for a given animal
#'   ID. Must include a column labeled \code{time1} that numbers each of the
#'   observations in consecutive order, which is automatically generated by
#'   \code{\link{filter_time}}.
#' @param brkpts A data frame of breakpoints for each animal ID (as generated by
#'   \code{\link{get_breakpts}}).
#'
#' @return A data frame that updates the original data object by including the
#'   segment number associated with each observation in relation to the extracted
#'   breakpoints.
#'
#'
#'
#'
assign_tseg_internal=function(dat, brkpts){

  tmp<- which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>%
    purrr::discard(is.na)
  breakpt<- as.numeric(breakpt[1,])

  breakpt1<- c(0, breakpt, Inf)
  n<- length(breakpt1)
  res<- matrix(NA, nrow(dat), 1)

  #in case time1 modified
  if (dat$time1[1] != 1 | dplyr::n_distinct(diff(dat$time1)) > 1) dat$time1<- seq(1, nrow(dat))

  for (i in 2:n){
    ind<- which(breakpt1[i-1] < dat$time1 & dat$time1 <= breakpt1[i])
    res[ind,]<- i-1
  }

  dat$tseg<- as.vector(res)
  dat
}
#------------------------------------------------

#' Add segment numbers to observations
#'
#' After breakpoints have been extracted for each animal ID, this function
#' assigns the associated segment number to observations for each animal ID.
#' These segments of observations will be used in the second stage of the model
#' framework to perform mixed-membership clustering by Latent Dirichlet
#' Allocation.
#'
#' @param dat A list where each element stores the data for a unique animal ID.
#'   Each element is a data frame that contains all data associated for a given
#'   animal ID and must include a column labeled \code{time1} that numbers each
#'   of the observations in consecutive order. This variable is automatically
#'   generated by the \code{\link{filter_time}} function during data
#'   preparation.
#' @param brkpts A data frame of breakpoints for each animal ID (as generated by
#'   \code{\link{get_breakpts}}).
#'
#' @return A data frame that updates the original data object by including the
#'   segment number associated with each observation in relation to the
#'   extracted breakpoints.
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
#' #future::plan(future::multisession)  #run all MCMC chains in parallel
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                            alpha = alpha)
#'
#'
#' # Determine MAP iteration for selecting breakpoints and store breakpoints
#' MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
#' brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)
#'
#'
#' # Assign track segments to all observations by ID
#' tracks.seg<- assign_tseg(dat = tracks.list, brkpts = brkpts)
#' }
#'
#' @export
assign_tseg=function(dat, brkpts){

  dat.out<- purrr::map(dat, assign_tseg_internal, brkpts = brkpts) %>%
    dplyr::bind_rows()

  dat.out
}
#------------------------------------------------

#' Convert data frame to a list by animal ID
#'
#' Converts an object of class \code{data.frame} to a list where each element is
#' a separate animal ID. This function prepares the data for further analysis
#' and when mapping other functions onto the data for separate animal IDs.
#'
#' @param dat A data frame containing the data for each animal ID.
#' @param ind character. The name of the column storing the animal IDs.
#'
#' @return A list where each element stores the data for a separate animal ID.
#'
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #convert to list
#' dat.list<- df_to_list(dat = tracks, ind = "id")
#'
#'
#' @export
df_to_list=function(dat, ind) {
  id<- dat %>%
    dplyr::pull(ind) %>%
    unique() %>%
    as.character()  #change to character in case of factors

  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id

  dat<- data.frame(dat)  #make sure as data.frame, not tibble
  dat[,ind]<- as.character(dat[,ind])  #change to character in case of factors

  for (i in 1:length(id)) {
    tmp<- which(dat[,ind] == id[i])
    dat.list[[i]]<- dat[tmp,]
  }
  dat.list
}
#------------------------------------------------

#' Find changes for integer variable
#'
#' Identify changes within a discrete variable. These values can be used to
#' pre-specify breakpoints within the segmentation model using
#' \code{\link{segment_behavior}}.
#'
#' @param dat A data frame containing the data for each animal ID.
#' @param ind character. The name of the column storing the discrete variable of
#'   interest.
#'
#' @return A vector of breakpoints is returned based on the data provided. If
#'   wishing to identify separate breakpoints per animal ID, this function
#'   should be mapped onto a list generated by \code{\link{df_to_list}}.
#'
#'
#' @examples
#' #simuluate data
#' var<- sample(1:3, size = 50, replace = TRUE)
#' var<- rep(var, each = 20)
#' id<- rep(1:10, each = 100)
#'
#' #create data frame
#' dat<- data.frame(id, var)
#'
#' #create list
#' dat.list<- df_to_list(dat = dat, ind = "id")
#'
#' #run function using purrr::map()
#' breaks<- purrr::map(dat.list, ~find_breaks(dat = ., ind = "var"))
#'
#' #or with lapply()
#' breaks1<- lapply(dat.list, find_breaks, ind = "var")
#'
#' @export
find_breaks=function(dat, ind) {
  which(diff(dat[,ind]) != 0) + 1
}
#------------------------------------------------

#' View trace-plots of output from Bayesian segmentation model
#'
#' Visualize trace-plots of the number of breakpoints estimated by the model as
#' well as the log marginal likelihood (LML) for each animal ID.
#'
#' @param data A list of model results that is returned as output from \code{\link{segment_behavior}}.
#' @param type character. The type of data that are being plotted from the
#'   Bayesian segmentation model results. Takes either 'nbrks' for the number of
#'   breakpoints or 'LML' for the log marginal likelihood.
#'
#' @return Trace-plots for the number of breakpoints or the log marginal
#'   likelihood are displayed for each of the animal IDs that were analyzed by
#'   the segmentation model.
#'
#'
#' @examples
#' \donttest{
#' #load data
#' data(tracks.list)
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
#' future::plan(future::multisession, workers = 3)  #run all MCMC chains in parallel
#'
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                                alpha = alpha)
#'
#' future::plan(future::sequential)  #return to single core
#'
#'
#' #run function
#' traceplot(data = dat.res, type = "nbrks")
#' traceplot(data = dat.res, type = "LML")
#' }
#'
#' @importFrom graphics "par"
#' @importFrom graphics "plot"
#' @importFrom graphics "abline"
#' @export
traceplot=function(data, type) {

  oldpar<- par(no.readonly = TRUE)  #store original parameters
  on.exit(par(oldpar))  #exit w/ original parameters

  identity<- names(data[["brkpts"]])

  ifelse(length(identity) == 1, par(mfrow=c(1,1)),
         ifelse(length(identity) == 2, par(mfrow=c(1,2)), par(mfrow=c(2,2))))

  for (i in 1:length(identity)) {
    par(ask=TRUE)
    ngibbs<- length(data[["brkpts"]][[i]])
    plot(x = 1:ngibbs,
         y = na.omit(as.numeric(data[[type]][i,-1])),
         type = "l",
         xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints",
                       ifelse(type == "LML","Log Marginal Likelihood",
                              stop("Need to select one of 'nbrks' or 'LML' for plotting"))),
         main = paste(identity[i]))
    abline(v = ngibbs/2, col = "red", lwd = 2)
  }
}

#---------------------------------------------

#' Internal function to find the maximum a posteriori (MAP) estimate of the MCMC
#' chain
#'
#' Internal function to be used by a wrapper.
#'
#' @param dat numeric. A vector of log marginal likelihood values for a given
#'   animal ID.
#' @param nburn numeric. The size of the burn-in phase after which the MAP
#'   estimate will be identified.
#'
#' @return A numeric value indicating the iteration after the burn-in phase that
#'   holds the MAP estimate.
#'
#'
#'
#' @export
get_MAP_internal=function(dat, nburn) {

  if (which.max(dat[-1]) < nburn) {
    MAP.est<- dat[-1] %>%
      order(decreasing = T)
    MAP.est<- MAP.est[MAP.est > nburn][1]
  } else {
    MAP.est<- which.max(dat[-1])
  }

  return(MAP.est)
}
#---------------------------------------------

#' Find the maximum a posteriori (MAP) estimate of the MCMC chain
#'
#' Identify the MCMC iteration that holds the MAP estimate. This will be used to
#' inform \code{\link{get_breakpts}} as to which breakpoints should be retained
#' on which to assign track segments to the observations of each animal ID.
#'
#' @param dat A data frame where each row holds the log marginal likelihood
#'   values at each iteration of the MCMC chain.
#' @param nburn numeric. The size of the burn-in phase after which the MAP
#'   estimate will be identified.
#'
#' @return A numeric vector of iterations at which the MAP estimate was found
#'   for each animal ID.
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
#' #future::plan(future::multisession)  #run all MCMC chains in parallel
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                            alpha = alpha)
#'
#'
#' # Determine MAP iteration for selecting breakpoints and store breakpoints
#' MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
#' }
#'
#' @export
get_MAP=function(dat, nburn) {
  MAP.est<- vector()
  for (i in 1:nrow(dat)) {
    MAP.est[i]<- get_MAP_internal(dat[i,], nburn)
  }

  MAP.est
}
#---------------------------------------------

#' Extract breakpoints for each animal ID
#'
#' @param dat A list of lists where animal IDs are separated as well as the
#'   breakpoints estimated for each iteration of the MCMC chain. This is stored
#'   within \code{breakpts} of model results returned after running
#'   \code{\link{segment_behavior}}.
#' @param MAP.est numeric. A vector of values at which the maximum a posteriori
#'   (MAP) estimate was identified for each of the animal IDs as returned by
#'   \code{\link{get_MAP}}. These must be in the same order as the data for the
#'   IDs supplied to \code{segment_behavior()}.
#'
#' @return A data frame where breakpoints are returned per animal ID within each
#'   row. For animal IDs that have fewer breakpoints than the maximum number
#'   that were estimated, \code{NA} values are used as place holders for these
#'   breakpoints that do not exist.
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
#' #future::plan(future::multisession)  #run all MCMC chains in parallel
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                            alpha = alpha)
#'
#'
#' # Determine MAP iteration for selecting breakpoints and store breakpoints
#' MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
#' brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)
#' }
#'
#' @export
get_breakpts=function(dat, MAP.est) {
  tmp<- list()

  for(i in 1:length(dat)) {
    ind<- MAP.est[i]
    tmp[[i]]<- dat[[i]][[ind]]
  }

  names(tmp)<- names(dat)
  max.length<- max(sapply(tmp, length))
  tmp<- lapply(tmp, function(x) { c(x, rep(NA, max.length-length(x)))}) %>%
    dplyr::bind_rows() %>%
    t() %>%
    data.frame()

  tmp<- cbind(id = names(dat), tmp)
  names(tmp)<- c('id', paste0("Brk_",1:(ncol(tmp)-1)))

  tmp
}
#------------------------------------------------

#' Internal function for plotting breakpoints over each of the data streams
#'
#' An internal function for plotting the results of the segmentation model.
#'
#' @param data A data frame for a single animal ID that contains columns for the
#'   ID, date or time variable, and each of the movement variables that were
#'   analyzed by \code{\link{segment_behavior}}. Data streams can be in
#'   continuous or discrete form.
#' @param as_date logical. If \code{TRUE}, plots breakpoints and data streams
#'   over the date. By default, this is set to \code{FALSE}.
#' @param var_names A vector of the column names for the movement variables to
#'   be plotted over time.
#' @param var_labels A vector of the labels to be plotted on the y-axis for each
#'   movement variable.  Set to \code{NULL} by default.
#' @param brkpts A data frame that contains the breakpoints associated with each
#'   animal ID. This data frame is returned by \code{\link{get_breakpts}}.
#'
#' @return A line plot for each movement variable showing how the estimated
#'   breakpoints relate to the underlying data. Depending on the user input for
#'   \code{var_names}, this may either be on the scale of the original
#'   continuous data or the discretized data.
#'
#'
#' @importFrom rlang .data
#'
plot_breakpoints_behav=function(data, as_date, var_names, var_labels, brkpts) {

  #index brkpts for particular id
  ind<- which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>%
    purrr::discard(is.na) %>%
    as.numeric()
  if (as_date == TRUE) breakpt<- data[breakpt, "date"]
  breakpt<- data.frame(breaks = breakpt)


  # Reformat data into long form using preferred time variable
  x<- ifelse(as_date == TRUE, "date", "time1")

  #in case time1 modified
  if (x == "time1" & dplyr::n_distinct(diff(data$time1)) > 1) data$time1<- seq(1, nrow(data))
  dat.long<- data[,c("id", x, var_names)] %>%
    tidyr::pivot_longer(cols = -c(1:2), names_to = "var", values_to = "value")

  # Store number of variables/data streams
  var.len<- length(var_names)

  #Relabel variables is var_labels specified
  if (!is.null(var_labels)) {
    for (i in 1:var.len) {
      dat.long$var<- gsub(x = dat.long$var, pattern = var_names[i], replacement = var_labels[i])
    }
  }



  print(
    ggplot(dat.long, aes(x=.data[[x]], y=.data$value, color=.data$var)) +
      geom_line(na.rm = TRUE, size = 0.25) +
      facet_wrap(~.data$var, scales = 'free', nrow = var.len, strip.position = "left") +
      scale_color_brewer("", palette = "Dark2", guide = FALSE) +
      geom_vline(data = breakpt, aes(xintercept = .data$breaks - 0.5),
                 color = "black", size = 0.7, alpha = 1) +
      labs(x = "\nTime", y = "") +
      ggtitle(paste(unique(data$id))) +
      theme_bw() +
      theme(axis.title = element_text(size = 18),
            axis.text = element_text(size = 12),
            strip.text = element_text(size = 14),
            strip.background = element_blank(),
            strip.placement = "outside",
            plot.title = element_text(size = 20),
            panel.grid.minor = element_blank())
  )
}
#------------------------------------------------

#' Plot breakpoints over a time series of each movement variable
#'
#' Visualize the breakpoints estimated by the segmentation model as they relate
#' to either the original (continuous) or discretized data. These plots assist
#' in determining whether too many or too few breakpoints were estimated as well
#' as whether the user needs to redefine how they discretized their data before
#' analysis.
#'
#' @param data A list where each element stores a data frame for a given animal
#'   ID. Each of these data frames contains columns for the ID, date or time1
#'   generated by \code{\link{filter_time}}, as well as each of the movement
#'   variables analyzed by \code{\link{segment_behavior}}.
#' @param as_date logical. If \code{TRUE}, plots breakpoints and data streams
#'   over the date. By default, this is set to \code{FALSE}.
#' @param var_names A vector of the column names for the movement variables to
#'   be plotted over time.
#' @param var_labels A vector of the labels to be plotted on the y-axis for each
#'   movement variable.  Set to \code{NULL} by default.
#' @param brkpts A data frame that contains the breakpoints associated with each
#'   animal ID. This data frame is returned by \code{\link{get_breakpts}}.
#'
#' @return A line plot per animal ID for each movement variable showing how the estimated
#'   breakpoints relate to the underlying data. Depending on the user input for
#'   \code{var_names}, this may either be on the scale of the original
#'   continuous data or the discretized data.
#'
#' @import ggplot2
#'
#' @importFrom graphics "par"
#' @examples
#'
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
#' #future::plan(future::multisession)  #run all MCMC chains in parallel
#' dat.res<- segment_behavior(data = tracks.list2, ngibbs = ngibbs, nbins = nbins,
#'                            alpha = alpha)
#'
#'
#' # Determine MAP iteration for selecting breakpoints and store breakpoints
#' MAP.est<- get_MAP(dat = dat.res$LML, nburn = ngibbs/2)
#' brkpts<- get_breakpts(dat = dat.res$brkpts, MAP.est = MAP.est)
#'
#'
#' #run function
#' plot_breakpoints(data = tracks.list, as_date = FALSE, var_names = c("step","angle"),
#'     var_labels = c("Step Length (m)", "Turning Angle (rad)"), brkpts = brkpts)
#' }
#'
#' @export
plot_breakpoints=function(data, as_date = FALSE, var_names, var_labels = NULL, brkpts) {

  oldpar<- par(no.readonly = TRUE)  #store original parameters
  on.exit(par(oldpar))  #exit w/ original parameters

  par(ask = TRUE)
  purrr::walk(data, ~plot_breakpoints_behav(., as_date = as_date, var_names = var_names,
                                            var_labels = var_labels, brkpts = brkpts))
}
#------------------------------------------------

#' Internal function to calculate step lengths, turning angles, and time steps
#'
#' An internal function that calculates step lengths, turning angles, and time
#' steps for a given animal ID.
#'
#' @param dat A data frame that contains the columns associated with the x and y
#'   coordinates as well as the date-time. For easier interpretation of the
#'   model results, it is recommended that coordinates be stored after UTM
#'   projection (meters) as opposed to unprojected in decimal degrees (map
#'   units). Date-time should be of class \code{POSIXct} and be labeled
#'   \code{date} within the data frame.
#' @param coord.names character. A vector of the column names under which the
#'   coordinates are stored. The name for the x coordinate should be listed
#'   first and the name for the y coordinate second.
#'
#' @return A data frame where all original data are returned and new columns are
#'   added for step length (\code{step}), turning angle (\code{angle}),
#'   net-squared displacement (\code{NSD}), and time step (\code{dt}).
#'
#'
#' @importFrom stats "na.omit"
#'
prep_data_internal=function(dat, coord.names) {

  #change names of coords to 'x' and 'y'
  dat<- dat %>%
    dplyr::rename(x = coord.names[1], y = coord.names[2])

  #calculate step length
  step<- sqrt((diff(dat[,"x"]))^2 + (diff(dat[,"y"]))^2)
  dat$step<- c(step, NA)

  #calculate turning angle
  angle<- vector()
  for (i in 2:nrow(dat)) {
    angle[i-1]<- atan2(dat[(i+1),"y"] - dat[i,"y"], dat[(i+1),"x"] - dat[i,"x"]) -
      atan2(dat[i,"y"] - dat[i-1,"y"], dat[i,"x"] - dat[i-1,"x"])
  }
  #adjust angle if |angle| > pi
  angle<- ifelse(angle > pi, angle - 2*pi,
                 ifelse(angle < -pi, 2*pi + angle, angle))
  dat$angle<- c(NA, angle)

  # calculate net-squared displacement
  x0<- dat[1,"x"]  #identify starting x-coord
  y0<- dat[1,"y"]  #identify starting y-coord
  displ<- sqrt((dat[,"x"] - x0)^2 + (dat[,"y"] - y0)^2)
  dat$NSD<- displ^2

  #calculate time steps
  dt<- difftime(dat$date, dplyr::lag(dat$date, 1), units = "secs") %>%
    na.omit() %>%
    as.numeric() %>%
    round()
  dat$dt<- c(dt, NA)



  dat
}
#------------------------------------------------

#' Calculate step lengths, turning angles, net-squared displacement, and time
#' steps
#'
#' Calculates step lengths, turning angles, and net-squared displacement based
#' on coordinates for each animal ID and calculates time steps based on the
#' date-time. Provides a self-contained method to calculate these variables
#' without needing to rely on other R packages (e.g., \code{adehabitatLT}).
#' However, functions from other packages can also be used to perform this step
#' in data preparation.
#'
#' @param dat A data frame that contains a column for animal IDs, the columns
#'   associated with the x and y coordinates, and a column for the date. For
#'   easier interpretation of the model results, it is recommended that
#'   coordinates be stored in a UTM projection (meters) as opposed to
#'   unprojected in decimal degrees (map units). Date-time should be of class
#'   \code{POSIXct} and be labeled \code{date} within the data frame.
#' @param coord.names character. A vector of the column names under which the
#'   coordinates are stored. The name for the x coordinate should be listed
#'   first and the name for the y coordinate second.
#' @param id character. The name of the column storing the animal IDs.
#'
#' @return A data frame where all original data are returned and new columns are
#'   added for step length (\code{step}), turning angle (\code{angle}),
#'   net-squared displacement (\code{NSD}), and time
#'   step (\code{dt}). Names for coordinates are changed to \code{x} and
#'   \code{y}. Units for \code{step} and \code{NSD} depend on the projection of the
#'   coordinates, \code{angle} is returned in radians, and \code{dt} is
#'   returned in seconds.
#'
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #subset only first track
#' tracks<- tracks[tracks$id == "id1",]
#'
#' #calculate step lengths and turning angles
#' tracks<- prep_data(dat = tracks, coord.names = c("x","y"), id = "id")
#'
#' @export
prep_data=function(dat, coord.names, id) {

  purrr::map(df_to_list(dat = dat, ind = id),
      ~prep_data_internal(., coord.names = coord.names)) %>%
    dplyr::bind_rows() %>%
    mutate_at(c("step","angle","NSD"), ~round(., 3))


}
#------------------------------------------------

#' Insert NA gaps to regularize a time series
#'
#' @param data A data frame that minimally contains columns for animal ID, date,
#'   and time step. These must be labeled \code{id}, \code{date}, and \code{dt},
#'   respectively, where date is of class \code{POSIXct}.
#' @param int integer. An integer that characterizes the desired interval on
#'   which to insert new rows.
#' @param units character. The units of the selected time interval \code{int},
#'   which can be selected from one of "secs", "mins", "hours", "days", or
#'   "weeks".
#'
#' @return A data frame where new rows have been inserted to regularize the \code{date} column. This results in values provided for \code{id}, \code{date}, and {dt} while inserting NAs for all other columns. Additionally, observations with duplicate date-times are removed.
#'
#' @examples
#' #load data
#' data(tracks)
#'
#' #remove rows to show how function works (create irregular time series)
#' set.seed(1)
#' ind<- sort(sample(2:15003, 500))
#'
#' tracks.red<- tracks[-ind,]
#'
#' #calculate step lengths, turning angles, net-squared displacement, and time steps
#' tracks.red<- prep_data(dat = tracks.red, coord.names = c("x","y"), id = "id")
#'
#' #round times to nearest interval
#' tracks.red<- round_track_time(dat = tracks.red, id = "id", int = c(3600, 7200, 10800, 14400),
#'                               tol = 300, units = "secs")
#'
#' #insert NA gaps
#' dat.out<- insert_NAs(tracks.red, int = 3600, units = "secs")
#'
#'
#' @export
insert_NAs = function(data, int, units) {

  dat.list<- bayesmove::df_to_list(data, "id")

  dat.list<- purrr::map(dat.list, ~{

    dat<- data.frame(.x)
    ind<- which(!is.na(dat$dt) & dat$dt > int)
    ind2<- ind + 1

    for (i in 1:length(ind)) {
      if (dat$dt[ind[i]] >= 2*int) {

        vec.length<- floor(dat$dt[ind[i]]/int)  #find multiple of int in dt to determine seq length
        seq.dates<- seq(dat$date[ind[i]], dat$date[ind2[i]],
                        by = paste(int, units))[1:vec.length]
        tmp1<- as.data.frame(lapply(dat[ind[i],], rep, length(seq.dates)))
        tmp1$date<- seq.dates
        tmp1$dt<- int
        NA.names<- which(!(names(tmp1) %in% c("id","date","dt")))
        tmp1[2:nrow(tmp1),NA.names]<- NA  #insert NAs for added obs

        dat[seq(ind[i]+vec.length, nrow(dat)+vec.length-1),]<- dat[ind2[i]:nrow(dat),]
        dat[ind[i]:(ind[i]+vec.length-1),]<- tmp1

        ind<- ind + (vec.length - 1)  #update index
        ind2<- ind2 + (vec.length - 1)  #update index

      } else {
        dat<- dat
      }
    }


    #update dt column and remove any duplicate rows
    dat$dt<- c(as.numeric(difftime(dat$date, dplyr::lag(dat$date), units = units))[-1], NA)
    dat<- dplyr::distinct(dat, date, .keep_all = TRUE)

    dat
  })

  dat.out<- dplyr::bind_rows(dat.list)

  dat.out
}
