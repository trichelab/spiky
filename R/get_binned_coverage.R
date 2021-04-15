#' tabulate read coverage in predefined bins
#'
#' refactored out of scan_spiked_bam
#'
#' @param bins    the GRanges with bins
#' @param covg    the coverage result (an RleList)
#'
#' @return        a GRanges of summarized coverage
#'
#' @seealso       get_spiked_coverage
#' @seealso       scan_spiked_bam
#'
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                    mustWork=TRUE)
#' si <- seqinfo_from_header(sb)
#' genome(si) <- "spike"
#' data(spike,package="spiky")
#' mgr <- get_merged_gr(si,spike)
#' fl <- scanBamFlag(isDuplicate=FALSE, isPaired=TRUE, isProperPair=TRUE)
#' bp <- ScanBamParam(flag=fl)
#' bamMapqFilter(bp) <- 20
#'
#' covg <- get_spiked_coverage(sb, bp=bp, gr=mgr)
#' get_binned_coverage(bins=GRanges(), covg=covg)
#'
#' @import        GenomicRanges
#'
#' @export
get_binned_coverage <- function(bins, covg) {

  if (length(bins) < 1) {
    message("Empty bins provided, skipping.")
    return(bins)
  }
  message("Binning genomic coverage...", appendLF=FALSE)
  bc <- binnedAverage(bins, numvar=covg[seqlevels(bins)], varname="coverage")
  message("Done.")
  return(bc)

}
