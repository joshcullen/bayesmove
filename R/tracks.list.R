#' Tracks discretized and prepared for segmentation.
#'
#' A dataset containing the prepared data after discretizing step lengths and
#' turning angles, as well as filtering observations at the primary time step.
#'
#' @format A list with three elements, each containing a data frame with ~4700
#'   rows and 11 variables: \describe{
#'   \item{id}{ID for each simulated track}
#'   \item{date}{date, recorded as datetime}
#'   \item{x}{x coordinate of tracks}
#'   \item{y}{y coordinate of tracks}
#'   \item{step}{the step length calculated as the distance between successive locations
#'   measured in units}
#'   \item{angle}{the relative turning angle measured in radians}
#'   \item{dt}{the time step or sampling interval between datetimes of successive
#'   observations}
#'   \item{obs}{the ordered number of observations per ID before filtering for the primary
#'   time step}
#'   \item{time1}{the ordered number of observations per ID after filtering for the primary
#'   time step}
#'   \item{SL}{discretized step lengths, separated into five bins}
#'   \item{TA}{discretized turning angles, separated into eight bins}
#'   }
#'
"tracks.list"
