#' get the (max, median, or mean) coverage for spike-in contigs from a BAM/CRAM
#'
#' @param covg      the coverage RleList
#' @param spike_gr  the spike-in GRanges 
#' @param spike     information about the spikes (default: load `spike`) 
#' @param how       how to summarize the per-spike coverage (max)
#' 
#' @return          a GRanges with summarized coverage and features for each
#' 
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky", 
#'                   mustWork=TRUE)
#' si <- seqinfo_from_header(sb) 
#' genome(si) <- "spike"
#' mgr <- get_merged_gr(si)
#'
#' fl <- scanBamFlag(isDuplicate=FALSE, isPaired=TRUE, isProperPair=TRUE)
#' bp <- ScanBamParam(flag=fl)
#' bamMapqFilter(bp) <- 20
#' 
#' data(spike, package="spiky") 
#' covg <- get_spiked_coverage(sb, bp=bp, gr=mgr)
#' get_spike_depth(covg, spike_gr=mgr, spike=spike)
#'
#' @export
get_spike_depth <- function(covg, spike_gr, spike, how=c("max", "mean")) {

  how <- match.fun(match.arg(how)) 
  if (!is(spike, "DFrame")) stop("Please provide a spike database")
  cols <- colnames(spike)[-1] 
  canon <- names(spike_gr)
  
  message("Summarizing spike-in counts...", appendLF=FALSE)
  spike_depth <- vapply(covg[seqlevels(spike_gr)], how, numeric(1), na.rm=TRUE)
  for (nm in cols) mcols(spike_gr)[[nm]] <- spike[canon, nm]
  spike_gr$coverage <- spike_depth[as.character(seqnames(spike_gr))]
  message("Done.") 
  
  return(spike_gr)

}
