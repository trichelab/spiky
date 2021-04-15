#' for spike-in contigs in GRanges, match to standardized spike seqlevels
#'
#' This function is essentially the opposite of rename_spikes, except that
#' it works well on GRanges/GAlignments from or for merged genome+spike BAMs.
#' If spike contigs are found, it will assign genome='spike' to those, while
#' changing the seqlevels to standardized names that match rownames(spike).
#'
#' @param x         something with seqlevels (GRanges, GAlignments, Seqinfo...)
#' @param spike     a DataFrame where spike$sequence is a DNAStringSet (or NULL)
#'
#' @return          x, but with standardized spike seqlevels and genomes
#'
#' @seealso         rename_spikes
#'
#' @import Rsamtools
#'
#' @export
rename_spike_seqlevels <- function(x, spike=NULL) {

  if (is.null(spike)) spike = spiky::spike
  has_spike <- find_spike_contigs(x, spike=spike)
  if (length(has_spike) > 0) {
    genome(x)[has_spike] <- "spike"
    new_contigs <- attr(has_spike, "mappings")
    old_contigs <- seqlevels(x)[has_spike]
    seqlevels(x)[has_spike] <- new_contigs[old_contigs]
    genome(x)[has_spike] <- "spike"
  } else {
    message("No spike-in contigs were found.")
  }
  return(x)

}
