[![Spiky](vignettes/spiky-small.png)](https://github.com/trichelab/spiky)

# spiky

[![Build Status](https://travis-ci.org/trichelab/spiky.png?branch=master)](https://travis-ci.org/trichelab/spiky)  [![codecov](https://codecov.io/gh/trichelab/spiky/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/spiky)

# Spike-ins for everyone 

Please note that `spiky` currently requires a recent version of R (4.0 or newer).
This restriction may be relaxed if we can verify it won't break things.

## Installing

The pre-release version of the package can be pulled from GitHub using the [BiocManager](https://cran.r-project.org/package/BiocManager) package:

    # install.packages("remotes")
    # install.packages("BiocManager")
    BiocManager::install("trichelab/spiky", build_vignettes=TRUE)

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
