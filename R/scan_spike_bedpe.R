#' Scan spikes BEDPE
#'
#'
#' @param bedpe       the BEDPE file path, or output from read_bedpe()
#' @param spike       information about the spikes (default: load `spike`)
#' @param how       how to summarize the per-spike coverage (max)
#' @return          a GRanges with coverage
#'
#' @examples
#' data(spike, package="spiky")
#' fl <- system.file("extdata", "example_spike_bedpe.bed.gz", package="spiky",mustWork=TRUE)
#' scan_spike_bedpe(fl,spike=spike) # will warn user about spike contigs
#'
#' @export
scan_spike_bedpe <- function(bedpe, spike, how="max"){

  if (!is(spike, "DFrame")) stop("Please provide a spike database")
  if (!is(bedpe, "Pairs")) {bedpe <- read_bedpe(bedpe)}
  s <- convertPairedGRtoGR(bedpe)
  seqlevels(s) <- sub("_mod_meth", "", seqlevels(s))
  seqlevels(s) <- sub("_meth", "", seqlevels(s))
  s_covg <- coverage(s)
  s_depth <- get_spike_depth(s_covg,spike=spike,how=how)
  names(s_depth) <- get_base_name(seqnames(s_depth))
  return(as(s_depth,"GRanges"))

}
