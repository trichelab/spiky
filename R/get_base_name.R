#' refactored out of rename_spikes and rename_spike_seqlevels
#'
#' A common task between generate_spike_fasta, rename_spikes, and
#' rename_spike_seqlevels is to determine what the largest common subset of
#' characters between existing contig names and stored standardized contigs
#' might be. This function eases that task.
#'
#' @param     contig_names    the names of contigs
#' @param     sep             separator character in contig names ("_")
#'
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' bh <- scanBamHeader(BamFile(sb))
#' orig_contigs <- names(bh$targets)
#' get_base_name(orig_contigs)
#'
#' @return                    a vector of elements 1:3 from each contig name
#'
#' @export
get_base_name <- function(contig_names, sep="_") {

  elts <- seq_len(3)
  contig_names <- as.character(contig_names)
  res <- vapply(lapply(strsplit(contig_names, sep), `[`, elts),
                paste, character(1), collapse=sep)
  names(res) <- contig_names
  #keep <- grep("NA", res, invert=TRUE)
  #res[keep]

}
