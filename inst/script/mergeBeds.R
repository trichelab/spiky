# This script turns a bunch of BED files into a GRanges with a matrix hanging off of it.

library(rtracklayer)

# in conjunction with tabixToGRanges, this can reduce arbitrarily many BEDs
mergeMcols <- function(x, y) {
  stopifnot (identical(granges(x), granges(y)))
  mcols(x)[, names(mcols(y))] <- mcols(y)
  return(x)
}

BEDs <- list.files(patt="bed.gz$")
names(BEDs) <- vapply(strsplit(BEDs, "\\."), `[`, 3, FUN.VALUE=character(1))

GRs <- lapply(BEDs, import)
for (samp in names(GRs)) names(mcols(GRs[[samp]]))[1] <- samp

merged <- Reduce(mergeMcols, GRs)
save(merged, file="../../Human01_pmol_test.rda")
