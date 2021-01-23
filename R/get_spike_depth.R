#' get the (max, median, or mean) coverage for spike-in contigs from a BAM/CRAM
#'
#' @param covg      the coverage RleList
#' @param spike_gr  the spike-in GRanges 
#' @param how       how to summarize the per-spike coverage (max)
#' @param spike     information about the spikes (default: load `spike`) 
#' 
#' @return          a GRanges with summarized coverage and features for each
#' 
#' @export
get_spike_depth <- function(covg, spike_gr, how=c("max", "mean"), spike=NULL) {

  how <- match.fun(match.arg(how)) 
  if (!is(spike, "DFrame")) data(spike)
  cols <- colnames(spike)[-1] 
  canon <- names(spike_gr)
  
  message("Summarizing spike-in counts...", appendLF=FALSE)
  spike_depth <- sapply(covg[seqlevels(spike_gr)], how, na.rm=TRUE)
  for (nm in cols) mcols(spike_gr)[[nm]] <- spike[canon, nm]
  spike_gr$coverage <- spike_depth[as.character(seqnames(spike_gr))]
  message("Done.") 
  
  return(spike_gr)

}
