[![Spiky](vignettes/spiky-small.png)](https://github.com/trichelab/spiky)

# spiky

[![Build Status](https://travis-ci.org/trichelab/spiky.png?branch=master)](https://travis-ci.org/trichelab/spiky)  [![codecov](https://codecov.io/gh/trichelab/spiky/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/spiky)

# Spike-ins for everyone 

This is a lightly edited README from skeletor which will be updated soon.

### Travis

Now you can go to [Travis](https://travis-ci.org/profile/trichelab) and turn on continuous integration for your new package. You may need to click the "Sync account" button to get your new package to show up in the list.

If you have a codecov.io account, running your tests on Travis will trigger the code coverage job. No additional configuration is necessary

## Installing

<!-- If you're putting `spiky` on CRAN, it can be installed with

    install.packages("spiky") -->

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
