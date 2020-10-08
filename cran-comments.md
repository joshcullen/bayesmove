## Resubmission
This is a resubmission. In this version I have:

* Added an on.exit() call immediately within functions where I adjust par() for plotting.

* Unwrapped some examples that used \dontrun{} and changed others to \donttest{} if they took > 5 s. However, the examples wrapped with \donttest{} were still checked by --run-donttest when pushed to Travis-CI and AppVeyor.

## Test environments
* local OS X install (and on travis-ci), R 4.0.2
* ubuntu-16.04 devel and release (on travis-ci and R-hub), R 4.0.2
* windows devel and release (on appveyor and win-builder), R 4.0.2

## R CMD check results
0 ERRORS | 0 WARNINGS | 1 NOTE

* checking checking CRAN incoming feasiblity ... NOTE
  Maintainer: 'Joshua Cullen <joshcullen10@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    biologging (14:71)
    pre (15:76)
  
  This is my first submission of the 'bayesmove' package to CRAN. Both 'biologging' and 'prep' are spelled correctly within the Description section of the DESCRIPTION file. The words 'biologging' and 'pre' are also saved within WORDLIST.
  

## Downstream dependencies
There are currently no downstream dependencies for this package.
