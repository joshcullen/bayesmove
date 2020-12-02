
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesmove <img src="man/figures/logo.png" align="right" width=35%  style="padding-left: 10px"/>

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/bayesmove)](https://CRAN.R-project.org/package=bayesmove)
[![Travis build
status](https://travis-ci.com/joshcullen/bayesmove.svg?branch=master)](https://travis-ci.com/joshcullen/bayesmove)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/joshcullen/bayesmove?branch=master&svg=true)](https://ci.appveyor.com/project/joshcullen/bayesmove)
[![Codecov test
coverage](https://codecov.io/gh/joshcullen/bayesmove/branch/master/graph/badge.svg)](https://codecov.io/gh/joshcullen/bayesmove?branch=master)
[![CRAN monthly
downloads](https://cranlogs.r-pkg.org/badges/bayesmove)](https://cran.r-project.org/web/packages/bayesmove/index.html)
[![CRAN total
downloads](https://cranlogs.r-pkg.org/badges/grand-total/bayesmove)](https://cran.r-project.org/web/packages/bayesmove/index.html)
<!-- badges: end -->

## Introduction

The goal of **bayesmove** is to analyze animal movement using a
non-parametric Bayesian framework, which addresses a number of
limitations of existing segmentation methods and state-space models.
This package currently offers two different model frameworks on which to
make behavioral inference from animal telemetry data: 1) **segment-level
behavioral state estimation** and 2) **observation-level behavioral
state estimation**.

The model that makes segment-level inference is a two-stage framework
that first partitions individual tracks into segments and subsequently
clusters these segments into latent behavioral states. This framework
allows the analysis of multiple telemetry and biologging data streams,
which must first be discretized into a set of bins before they can be
analyzed. The model that makes behavioral inference at the
observation-level also requires that data streams are first discretized,
but then directly clusters these observations together into behavioral
states within a single step. While the outcome is similar to that from
state-space and hidden Markov models, this observation-level model does
not assume an underlying Markov property or use a mechanistic process
(e.g., correlated random walk).

This package also includes features to check model convergence based on
the log-likelihood for each MCMC iteration. Model output are often
returned in a format that is `tidyverse`-friendly, which allows for easy
visualization using `ggplot2`.

## Installation

You can install the latest CRAN release with:

``` r
install.packages("bayesmove")
```

You can install the latest **stable** version of the package from GitHub
with:

``` r
# install.packages("remotes")
remotes::install_github("joshcullen/bayesmove")
```

or latest development (**unstable**) version with:

``` r
# install.packages("remotes")
remotes::install_github("joshcullen/bayesmove@dev")
```

## Support

If you are receiving errors from the model output that you believe to be
bugs, please report them as issues in the GitHub repo. Additionally, if
there are any other features you would like added to this package,
please submit them to the issue tracker.
