% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/segmentation_helper_functions.R
\name{plot_breakpoints_behav}
\alias{plot_breakpoints_behav}
\title{Internal function for plotting breakpoints over each of the data streams}
\usage{
plot_breakpoints_behav(data, as_date, var_names, var_labels, brkpts)
}
\arguments{
\item{data}{A data frame for a single animal ID that contains columns for the
ID, date or time variable, and each of the movement variables that were
analyzed by \code{\link{segment_behavior}}. Data streams can be in
continuous or discrete form.}

\item{as_date}{logical. If \code{TRUE}, plots breakpoints and data streams
over the date. By default, this is set to \code{FALSE}.}

\item{var_names}{A vector of the column names for the movement variables to
be plotted over time.}

\item{var_labels}{A vector of the labels to be plotted on the y-axis for each
movement variable.  Set to \code{NULL} by default.}

\item{brkpts}{A data frame that contains the breakpoints associated with each
animal ID. This data frame is returned by \code{\link{get_breakpts}}.}
}
\value{
A line plot for each movement variable showing how the estimated
  breakpoints relate to the underlying data. Depending on the user input for
  \code{var_names}, this may either be on the scale of the original
  continuous data or the discretized data.
}
\description{
An internal function for plotting the results of the segmentation model.
}
