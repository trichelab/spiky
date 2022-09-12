#' compute methylation specificity for spike-in standards
#'
#' In a cfMeDIP experiment, the yield of methylated fragments should be >95%
#' (ideally 98-99%) due to the nature of the assay.
#'
#' @param   spike_gr    GRanges of spike contigs (e.g. output object from scan_spiked_bam, scan_spike_contigs, or scan_spike_bedpe)
#' @param   spike       spike contig database, if needed (e.g. data(spike))
#'
#' @return              list with median and mean coverage across spike contigs
#'
#' @examples
#' data(genomic_res)
#' data(spike_res)
#' data(spike, package="spiky")
#' methylation_specificity(spike_res, spike=spike)
#'
#' @export
methylation_specificity <- function(spike_gr, spike) {
  methreads <- covg_to_df(spike_gr, meth=TRUE, spike=spike)$read_count
  totalreads <- covg_to_df(spike_gr, meth=FALSE, spike=spike)$read_count
  meth_spec <- list("mean" = base::mean(methreads/totalreads),
                    "median" = stats::median(methreads/totalreads))
  return(meth_spec)
}
