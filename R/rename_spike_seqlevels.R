#' for spike-in contigs in GRanges, match to standardized spike seqlevels
#' 
#' This function is essentially the opposite of rename_spikes, except that 
#' it works well on GRanges/GAlignments from or for merged genome+spike BAMs.
#' 
#' @param x         something with seqlevels
#' @param spike     a DataFrame where spike$sequence is a DNAStringSet
#' 
#' @return          x, but with standardized spike seqlevels
#' 
#' @seealso         rename_spikes
#' 
#' @import Rsamtools
#' @import Biostrings
#' 
#' @export
rename_spike_seqlevels <- function(x, spike) { 

  orig_contigs <- seqlevels(x)
  new_contigs <- get_base_name(seqlevels(x))
  has_spike <- which(new_contigs %in% rownames(spike))
  seqlevels(x)[has_spike] <- new_contigs[has_spike] 
  return(x)

}
