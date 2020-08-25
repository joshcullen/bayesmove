

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
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' id<- rep(1:10, each = 100)
#'
#' #create data frame
#' dat<- data.frame(id, step, angle)
#'
#' #define limits for each bin
#' dist.lims<- c(quantile(step, c(0, 0.25, 0.5, 0.75, 0.95)), max(step))
#' angle.lims<- c(-pi, -3*pi/4, -pi/2, -pi/4, 0, pi/4, pi/2, 3*pi/4, pi)
#'
#' #run function
#' dat1<- discrete_move_var(dat = dat, lims = list(dist.lims, angle.lims),
#'                          varIn = c("step", "angle"),
#'                          varOut = c("SL","TA"))
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
#'
#' @return A data frame where \code{dt} and \code{date} are both adjusted based
#'   upon the rounding of time intervals according to the specified tolerance.
#'
#'
#' @examples
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
#' date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
#' dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
#' dt<- c(dt, NA)
#' id<- rep(1:10, each = 100)
#'
#' #create data frame
#' dat<- data.frame(id, date, dt, step, angle)
#'
#' #run function
#' dat1<- round_track_time(dat = dat, id = "id", int = 3600, tol = 20, time.zone = "UTC")
#'
#' @export
round_track_time = function(dat, id, int, tol, time.zone = "UTC") {

  dat<- df_to_list(dat, ind = id)
  for (i in 1:length(dat)) {
    tmp<- matrix(NA, nrow(dat[[i]]), 2)

    if (length(int) == 1) {  #when using only 1 time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j, 1:2]<- NA
        } else if (dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                   dat[[i]]$dt[j] != int) {
          tmp[j, 1:2]<- c(int, dat[[i]]$date[j] - (dat[[i]]$dt[j] - int))
        } else {
          tmp[j, 1:2]<- c(dat[[i]]$dt[j], dat[[i]]$date[j])
        }
      }
    } else {  #when using more than one time interval
      for (j in 1:nrow(dat[[i]])) {
        if (is.na(dat[[i]]$dt[j])) {
          tmp[j, 1:2]<- NA
        } else if (sum(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                       dat[[i]]$dt[j] != int) > 0) {
          ind<- which(dat[[i]]$dt[j] >= (int - tol) & dat[[i]]$dt[j] <= (int + tol) &
                        dat[[i]]$dt[j] != int)
          tmp[j, 1:2]<- c(int[ind], dat[[i]]$date[j] - (dat[[i]]$dt[j] - int))
        } else {
          tmp[j, 1:2]<- c(dat[[i]]$dt[j], dat[[i]]$date[j])
        }
      }
    }
    dat[[i]]$dt<- tmp[,1]
    dat[[i]]$date<- tmp[,2] %>%
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
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
#' date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
#' dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
#' dt<- c(dt, NA)
#' id<- rep(1:10, each = 100)
#'
#' #create data frame
#' dat<- data.frame(id, date, dt, step, angle)
#' dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC")
#'
#' #create list
#' dat.list<- df_to_list(dat = dat, ind = "id")
#'
#' #run function
#' dat.list.filt<- filter_time(dat.list = dat.list, int = 3600)
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
    behav.list[[i]]<- dat.list[[i]] %>%
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
#' @export
behav_seg_image=function(dat, nbins) {

  dat<- dat[,-1]  #remove id col
  behav.list<- purrr::map2(list(dat), nbins, ~matrix(0, nrow = nrow(.x), ncol = .y))
  for (i in 1:length(behav.list)) {
    for (j in 1:nrow(dat)){
      behav.list[[i]][j,dat[,i][j]]=1
    }
  }

  names(behav.list)<- names(dat)
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
#' @export
assign_tseg_internal=function(dat, brkpts){

  tmp<- which(unique(dat$id) == brkpts$id)
  breakpt<- brkpts[tmp,-1] %>%
    purrr::discard(is.na)
  breakpt<- as.numeric(breakpt[1,])

  breakpt1<- c(0, breakpt, Inf)
  n<- length(breakpt1)
  res<- matrix(NA, nrow(dat), 1)
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
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
#' date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
#' dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
#' dt<- c(dt, NA)
#' id<- rep(1:10, each = 100)
#'
#' #simulate breakpoints
#' brkpts<- rep(sort(sample(1:65, 7, replace = TRUE)), 10)
#' brkpts<- matrix(brkpts, 10, 7, byrow = TRUE)
#' brkpts<- data.frame(id = 1:10, brkpts)
#'
#' #create data frame
#' dat<- data.frame(id, date, dt, step, angle)
#' dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC")
#'
#' #create list
#' dat.list<- df_to_list(dat = dat, ind = "id")
#'
#' #filter by primary time step
#' dat.list.filt<- filter_time(dat.list = dat.list, int = 3600)
#'
#' #run function
#' dat1<- assign_tseg(dat = dat.list.filt, brkpts = brkpts)
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
#' #simulate data
#' step<- rgamma(1000, c(1, 2.5, 10), c(1, 1, 1))
#' angle<- runif(1000, -pi, pi)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 1000)
#' date<- date + lubridate::seconds(runif(length(date), -15, 15))  #introduce noise
#' dt<- as.numeric(diff(date)) * 60  #convert time difference to seconds
#' dt<- c(dt, NA)
#' id<- rep(1:10, each = 100)
#'
#' #create data frame
#' dat<- data.frame(id, date, dt, step, angle)
#' dat<- round_track_time(dat = dat, id = "id", int = 3600, tol = 15, time.zone = "UTC")
#'
#' #run function
#' dat.list<- df_to_list(dat = dat, ind = "id")
#'
#'
#' @export
df_to_list=function(dat, ind) {
  id<- dat %>%
    dplyr::pull(ind) %>%
    unique()

  n=length(id)
  dat.list<- vector("list", n)
  names(dat.list)<- id

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
#' @param data A data frame containing values at each iteration of the MCMC
#'   chain where each row holds data for a specific animal ID.
#' @param ngibbs numeric. The total number of iterations of the MCMC chain.
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
#' #simulate data
#' ngibbs<- 1000
#' y<- (-1000 * 501:1500)/(-500 + 501:1500) + rnorm(ngibbs, 0, 0.1)
#' dat<- matrix(c(1, y), 1, 1001)
#' dat<- data.frame(dat)
#' names(dat)[1]<- "id"
#'
#' #run function
#' traceplot(data = dat, ngibbs = ngibbs, type = "LML")
#'
#' @importFrom graphics "par"
#' @importFrom graphics "plot"
#' @export
traceplot=function(data, ngibbs, type) {

  identity<- data$id

  ifelse(length(identity) == 1, par(mfrow=c(1,1)),
         ifelse(length(identity) == 2, par(mfrow=c(1,2)), par(mfrow=c(2,2))))

  for (i in 1:length(identity)) {
    par(ask=TRUE)
    plot(x=1:ngibbs, y=data[i,-1], type = "l",
         xlab = "Iteration",
         ylab = ifelse(type == "nbrks", "# of Breakpoints",
                       ifelse(type == "LML","Log Marginal Likelihood",
                              stop("Need to select one of 'nbrks' or 'LML' for plotting"))),
         main = paste(identity[i]))
  }
  on.exit(par(ask = FALSE, mfrow=c(1,1)))
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
#' #simulate data
#' ngibbs<- 1000
#' y<- (-1000 * 501:1500)/(-500 + 501:1500) + rnorm(ngibbs, 0, 0.1)
#' dat<- matrix(c(1, y), 1, 1001)
#' dat<- data.frame(dat)
#' names(dat)[1]<- "id"
#'
#' #run function
#' MAP.est<- get_MAP(dat = dat, nburn = ngibbs/2)
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
#' #simulate data
#' id1 = id2 = id3 = list(4, c(2,4), 4, 4, c(4,8), c(4,8,17), c(4,8,17),
#'                        c(4,8,20), c(4,8,20,25), c(5,8,20,25))
#' dat.list<- list(id1 = id1, id2 = id2, id3 = id3)
#'
#' MAP.est<- c(5, 8, 9)
#'
#'
#' #run function
#' dat1<- get_breakpts(dat = dat.list, MAP.est = MAP.est)
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

#' Internal function for plotting breakpoints over a heatmap of the discretized
#' movement variables
#'
#' An internal function for plotting the results of the segmentation model.
#'
#' @param data A data frame for a single animal ID that contains only columns
#'   for the ID and each of the movement variables that were analyzed by
#'   \code{\link{segment_behavior}}. The ID column must be first.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{data}.
#' @param brkpts A data frame that contains the breakpoints associated with each
#'   animal ID. This data frame is returned by \code{\link{get_breakpts}}.
#' @param title logical. If \code{TRUE}, includes the animal ID as the title of
#'   the plot.
#' @param legend logical. If \code{TRUE}, shows the legend for the heatmap.
#'
#' @return A heatmap for each movement variable showing the bins assigned to
#'   each observation over time that is overlayed by the extracted breakpoints.
#'
#'
#' @importFrom rlang .data
#' @export
plot_heatmap_behav=function(data, nbins, brkpts, title, legend) {

  #transform into pres/abs matrix
  behav.heat<- behav_seg_image(data, nbins)

  #convert to long form
  behav.heat.long<- purrr::map2(behav.heat, as.list(names(behav.heat)), ~{
    tmp<- .x %>%
      data.frame()
    tmp %>%
      tidyr::pivot_longer(cols = 1:ncol(tmp), names_to = "bin", values_to = "value") %>%
      dplyr::mutate(time = rep(1:nrow(data), each = length(unique(.data$bin))),
                    var = rep(.y, length(.data$bin)))}
  ) %>%
    dplyr::bind_rows() %>%
    dplyr::mutate_at("value", factor) %>%
    dplyr::mutate_at("bin", readr::parse_number)

  levels(behav.heat.long$value)<- c("Unused","Used")

  #index brkpts for particular id
  ind<- which(unique(data$id) == brkpts$id)
  breakpt<- brkpts[ind,-1] %>%
    purrr::discard(is.na) %>%
    t() %>%
    data.frame()
  names(breakpt)<- "breaks"


  if(legend == TRUE) {
    legend.pos<- "top"
  } else {
    legend.pos<- "none"
  }
  title<- ifelse(title == TRUE,
                 list(theme(axis.title = element_text(size = 18),
                            axis.text = element_text(size = 12),
                            strip.text = element_text(size = 12, face = 'bold'),
                            plot.title = element_text(size = 20, hjust = 0, vjust = -6),
                            plot.margin = margin(0, 1, 0.5, 0.5, "cm"),
                            legend.justification = "right",
                            legend.position = legend.pos,
                            legend.text = element_text(
                              margin = margin(r = 15, unit = "pt")))),
                 list(theme(axis.title = element_text(size = 18),
                            axis.text = element_text(size = 12),
                            strip.text = element_text(size = 12, face = 'bold'),
                            plot.title = element_blank(),
                            legend.justification = "right",
                            legend.position = legend.pos,
                            legend.text = element_text(
                              margin = margin(r = 15, unit = "pt")))))

  print(
    ggplot(behav.heat.long, aes(x=.data$time, y=as.character(.data$bin), fill=.data$value)) +
      geom_tile() +
      facet_wrap(~.data$var, scales = 'free', nrow = 2) +
      scale_fill_viridis_d('') +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      geom_vline(data = breakpt, aes(xintercept = .data$breaks - 0.5),
                 color = viridis::viridis(n=9)[7], size = 0.6, alpha = 1) +
      labs(x = "\nTime", y = "Bin\n") +
      ggtitle(paste(unique(data$id))) +
      theme_bw() +
      title
  )
}
#------------------------------------------------

#' Plot breakpoints over a time series heatmap of the movement variables
#'
#' Visualize the breakpoints estimated by the segmentation model as they relate
#' to the discretized data on which the model was performed. These plots assist
#' in determining whether too many or too few breakpoints were estimated as well
#' as whether the user needs to redefine how they discretized their data before
#' analysis.
#'
#' @param data A list where each element stores a data frame for a given animal
#'   ID. Each of these data frames contain only columns for the ID and each of
#'   the movement variables that were analyzed by
#'   \code{\link{segment_behavior}}. The ID column must be first.
#' @param nbins numeric. A vector of the number of bins used to discretize each
#'   movement variable. These must be in the same order as the columns within
#'   \code{data}.
#' @param brkpts A data frame that contains the breakpoints associated with each
#'   animal ID. This data frame is returned by \code{\link{get_breakpts}}.
#' @param title logical. If \code{TRUE}, includes the animal ID as the title of
#'   the plot.
#' @param legend logical. If \code{TRUE}, shows the legend for the heatmap.
#'
#' @return A heatmap per animal ID that shows the value of each movement
#'   variable over time (as discretized into bins) and is overlayed by the
#'   extracted breakpoints.
#'
#' @import ggplot2
#'
#' @importFrom graphics "par"
#' @examples
#'
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
#' dat.list.filt1<- lapply(dat.list.filt,
#'                         function(x) subset(x, select = c(id, SL, TA)))
#'
#'
#' #simulate breakpoints
#' brkpts<- rep(sort(sample(1:65, 7, replace = TRUE)), 10)
#' brkpts<- matrix(brkpts, 10, 7, byrow = TRUE)
#' brkpts<- data.frame(id = 1:10, brkpts)
#'
#'
#' #run function
#' plot_heatmap(data = dat.list.filt1, nbins = c(5,8), brkpts = brkpts,
#'              title = TRUE, legend = TRUE)
#'}
#'
#' @export
plot_heatmap=function(data, nbins, brkpts, title, legend) {

  par(ask = TRUE)
  purrr::map(data, ~plot_heatmap_behav(., nbins = nbins, brkpts = brkpts, title = title,
                                legend = legend))
  par(ask = FALSE)

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
#'   added for step length (\code{step}), turning angle (\code{angle}), and time
#'   step (\code{dt}).
#'
#'
#' @importFrom stats "na.omit"
#' @export
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

  #calculate time steps
  dt<- difftime(dat$date, dplyr::lag(dat$date, 1), units = "secs") %>%
    na.omit() %>%
    as.numeric() %>%
    round()
  dat$dt<- c(dt, NA)

  dat
}
#------------------------------------------------

#' Calculate step lengths, turning angles, and time steps
#'
#' Calculates step lengths and turning angles based on coordinates for each
#' animal ID and calculates time steps based on the date-time. Provides a
#' self-contained method to calculate these variables without needing to rely on
#' other R packages (e.g., \code{adehabitatLT}). However, functions from other
#' packages can also be used to perform this step in data preparation.
#'
#' @param dat A data frame that contains a column for animal IDs, the columns
#'   associated with the x and y coordinates, and a column for the date. For
#'   easier interpretation of the model results, it is recommended that
#'   coordinates be stored after UTM projection (meters) as opposed to
#'   unprojected in decimal degrees (map units). Date-time should be of class
#'   \code{POSIXct} and be labeled \code{date} within the data frame.
#' @param coord.names character. A vector of the column names under which the
#'   coordinates are stored. The name for the x coordinate should be listed
#'   first and the name for the y coordinate second.
#' @param id character. The name of the column storing the animal IDs.
#'
#' @return A data frame where all original data are returned and new columns are
#'   added for step length (\code{step}), turning angle (\code{angle}), and time
#'   step (\code{dt}). Names for coordinates are changed to \code{x} and
#'   \code{y}. Units for step length depend on the projection of the
#'   coordinates, turning angles are returned in radians, and time steps are
#'   returned in seconds.
#'
#'
#' @examples
#' #simulate data
#' lon<- c(1,1,3,4,4,5,7,9,10,13)
#' lat<- c(2,1,2,2,3,5,8,1,1,2)
#' date<- seq(c(ISOdate(2020, 6, 17, tz = "UTC")), by = "hour", length.out = 10)
#' date<- date + lubridate::seconds(runif(length(date), -300, 300))  #introduce noise
#' dat<- data.frame(id = 1, date, lon, lat)
#'
#' #run function
#' dat1<- prep_data(dat = dat, coord.names = c("lon","lat"), id = "id")
#'
#' @export
prep_data=function(dat, coord.names, id) {

  purrr::map(df_to_list(dat = dat, ind = id),
      ~prep_data_internal(., coord.names = coord.names)) %>%
    dplyr::bind_rows() %>%
    mutate_at(c("step","angle"), ~round(., 3))


}
