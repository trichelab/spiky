#' find spike-in seqlevels in an object `x`, where !is.null(seqinfo(x))
#' 
#' Find the spike-like contigs in a BAM with both natural and spiked contigs.
#' This started out as glue in some other functions and got refactored out. 
#' 
#' @param x       something with seqlevels
#' @param spike   a DataFrame with spike-in information, if not using defaults
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
#' find_spike_contigs(si)
#' 
#' @seealso get_base_name
#' @seealso rename_spike_seqlevels
#' 
#' @export
find_spike_contigs <- function(x, spike=NULL) { 

  if (is.null(spike)) {
    rm(spike)
    data(spike, package="spiky")
  }

  orig_contigs <- seqlevels(x)
  new_contigs <- get_base_name(seqlevels(x))
  res <- which(new_contigs %in% rownames(spike))

  mappings <- new_contigs
  names(mappings) <- orig_contigs
  attr(res, "mappings") <- mappings[res]
  return(res) 

}
