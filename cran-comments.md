## Resubmission
This is a resubmission. In this version I have:

* Fixed a note related to a possibly invalid file URI (URI: reference/figures/logo.png) from the README.md file by changing the image link to man/figures/logo.png. This NOTE disappeared when running again on win-builder and travis-ci.

## Test environments
* local OS X install (and on travis-ci), R 4.0.2
* ubuntu-16.04 devel and release (on travis-ci and R-hub), R 4.0.2
* windows (on appveyor and win-builder), R 4.0.2

## R CMD check results
0 ERRORS | 0 WARNINGS | 2 NOTES

* checking checking CRAN incoming feasiblity ... NOTE
  Maintainer: 'Joshua Cullen <joshcullen10@gmail.com>'

  New submission

  Possibly mis-spelled words in DESCRIPTION:
    biologging (14:71)
    pre (15:76)
  
  This is my first submission of the 'bayesmove' package to CRAN. Both 'biologging' and 'prep' are spelled correctly within the Description section of the DESCRIPTION file. The words 'biologging' and 'pre' are also saved within WORDLIST.
  

* checking for future file timestamps ... NOTE
  unable to verify current time
  
  For this note, the external clock used for this check appears to be down (<https://stackoverflow.com/questions/63613301/r-cmd-check-note-unable-to-verify-current-time>)

## Downstream dependencies
There are currently no downstream dependencies for this package.