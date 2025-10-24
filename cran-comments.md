## Test environments
* local macOS 15.3.2, R 4.4.2
* macos-15-arm64, R 4.5.1 on GitHub Actions
* ubuntu-24.04, R 4.4.3, 4.5.1, devel on GitHub Actions
* windows-2025, R 4.5.1 on GitHub Actions

## R CMD check results
0 ERRORS | 0 WARNINGS | 0 NOTES

Fixed outstanding NOTES listed on CRAN. This is my fourth submission of the 'bayesmove' package to CRAN.

## Summary of updates in v0.2.3

* Greatly improved speed and efficiency of get_MAP() function.
* Updated shiny_tracks() function by adding option to color mapped points by a selected variable in a dropdown menu.
* Fixed issue with expand_behavior() function where wrong sample size was specified for track segments by ID.
* Updated suggested citations for use of package.
* Fixed bug with use of assign_behavior().

## Revdep check results
There are currently no reverse dependencies for this package.
