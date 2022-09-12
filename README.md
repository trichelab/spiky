[![Spiky](vignettes/spiky-small.png)](https://github.com/trichelab/spiky)

# spiky

[![Build Status](https://travis-ci.org/trichelab/spiky.png?branch=master)](https://travis-ci.org/trichelab/spiky)  [![codecov](https://codecov.io/gh/trichelab/spiky/branch/master/graph/badge.svg)](https://codecov.io/gh/trichelab/spiky)

# Spike-ins for everyone 
Analysis of cfMeDIP-seq data with spike-in control sequences. For more information, see Wilson, et. al, Sensitive and reproducible cell-free methylome quantification with synthetic spike-in controls, Cell Reports Methods, 2022. (https://www.sciencedirect.com/science/article/pii/S266723752200176X?via%3Dihub)

Please note that `spiky` strongly suggests a recent version of R (4.0 or newer; 3.6 or newer is required).

This restriction may be removed if our dependencies can tolerate it.

## Installing

To install this package, start R (version "3.6" or later) and enter:

```
if (!requireNamespace("BiocManager", quietly = TRUE)) { 
  install.packages("BiocManager")
}
BiocManager::install("trichelab/spiky")
```

## For developers

The repository includes a Makefile to facilitate some common tasks.

### Running tests

`$ make test`. Requires the [testthat](https://github.com/hadley/testthat) package. You can also specify a specific test file or files to run by adding a "file=" argument, like `$ make test file=logging`. `test_package` will do a regular-expression pattern match within the file names. See its documentation in the `testthat` package.

### Updating documentation

`$ make doc`. Requires the [roxygen2](https://github.com/klutometis/roxygen) package.
