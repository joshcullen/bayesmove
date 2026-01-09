## Test environments
* local macOS 15.3.2, R 4.5.2
* macos-15-arm64, R 4.5.2 on GitHub Actions
* ubuntu-24.04, R 4.4.3, 4.5.2, devel on GitHub Actions
* windows-2025, R 4.5.2 on GitHub Actions

## R CMD check results
0 ERRORS | 0 WARNINGS | 0 NOTES

This package was flagged during a reverse dependency check by maintainer of {dplyr} due to a function that will no longer be exported. An update was made to {bayesmove} to accommodate this change by {dplyr}. This is my fourth submission of the {bayesmove} package to CRAN.

## Summary of updates in v0.2.4

* Now explicitly sets "id" as global variable to mitigate breaking change by {dplyr}.

## Revdep check results
There are currently no reverse dependencies for this package.
