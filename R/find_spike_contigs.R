#' find spike-in seqlevels in an object `x`, where !is.null(seqinfo(x))
#'
#' Find the spike-like contigs in a BAM with both natural and spiked contigs.
#' This started out as glue in some other functions and got refactored out.
#'
#' @param x       something with seqlevels
#' @param spike   a DataFrame with spike-in information
#'
#' @return        indices of which contigs in seqlevels(x) are spike-in contigs
#'
#' @details
#' The indices have an attribute "mappings", which is a character vector
#' such that attr(find_spike_contigs(x), "mappings") == standardized
#' for all contig names in the CRAM/BAM/whatever, and standardized is
#' the rowname in `spike` that corresponds to the original contig name.
#'
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' si <- seqinfo_from_header(sb)
#' data(spike, package="spiky")
#' find_spike_contigs(si, spike=spike)
#'
#' @seealso get_base_name
#' @seealso rename_spike_seqlevels
#'
#' @export
find_spike_contigs <- function(x, spike) {

  orig_contigs <- seqlevels(x)
  new_contigs <- get_base_name(seqlevels(x))
  #new_contigs <- orig_contigs
  res <- which(new_contigs %in% rownames(spike))

  mappings <- new_contigs
  if (length(new_contigs) == length(orig_contigs)) {
    names(mappings) <- orig_contigs
  } else {
    warning("Cannot resolve contig name mismatches with spike database names. Your file may not have any spikes, or the spike database might not contain them.\n")
  }

  attr(res, "mappings") <- mappings[res]
  return(res)

}
