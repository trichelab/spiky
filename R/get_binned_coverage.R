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
#' @import        GenomicRanges
#'
#' @export 
get_binned_coverage <- function(bins, covg) { 

  message("Binning genomic coverage...", appendLF=FALSE)
  bc <- binnedAverage(bins, numvar=covg[seqlevels(bins)], varname="coverage")
  message("Done.") 
  return(bc) 

}
