## Test environments
* local OS X install (and using github actions), R 4.1.1
* ubuntu-20.04 devel and release (using github actions and R-hub), R 4.1.1
* windows devel and release (using github actions and win-builder), R 4.1.1

## R CMD check results
0 ERRORS | 0 WARNINGS | 1 NOTE

* Possibly misspelled words in DESCRIPTION:
  al (19:64)
  et (19:61)

The words listed here are spelled correctly and are also included in the WORDLIST file. This is my third submission of the 'bayesmove' package to CRAN.

## Summary of updates in v0.2.1

* Greatly improved speed of segmentation model by adding internal function summarize1 in C++ and updated get_summary_stats().
* Updated shiny_tracks by improving responsiveness to changes in the time window, adding a data table with options to filter, and the ability to explore the time series and map of multiple IDs at once.
* Updated traceplot() to automatically determine the number of MCMC iterations and add a line denoting the burn-in period.
* Updated progress bar in segment_behavior() to prevent it from disappearing.

## Downstream dependencies
There are currently no downstream dependencies for this package.
