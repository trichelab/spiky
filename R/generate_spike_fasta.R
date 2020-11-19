#' for CRAM files, a FASTA reference is required to decode; this builds that
#' 
#' If the contigs in a CRAM have even slightly different names from those in
#' the reference, decoding will fail.  In some cases there are multiple names
#' for a given contig (which raises the question of whether to condense them),
#' and thus the same reference sequence decodes multiple contig names. 
#' 
#' This function simply generates an appropriate spike reference for a CRAM,
#' using the CRAM headers to figure out which references are used for which. 
#' At the moment, CRAM support in Rsamtools only exists in the GitHub branch:
#' 
#' BiocManager::install("Bioconductor/Rsamtools@cram") 
#'
#' Using other versions of Rsamtools will yield an error. 
#' 
#' @param x         a CRAM file, hopefully with an index
#' @param fasta     the filename for the resulting FASTA ("spikes.fa")
#' @param spike     a DataFrame where spike$sequence is a DNAStringSet (spike)
#' 
#' @return          a DNAStringSet with renamed contigs, as exported to `fasta` 
#' 
#' @import Biostrings 
#' @import Rsamtools 
#' 
#' @export
generate_spike_fasta <- function(x, fasta="spike_contigs.fa", spike=NULL) { 

  cram <- BamFile(x) 
  hdr <- scanBamHeader(cram) 
  if (is.null(spike)) data(spike, package="spiky") 
  cram_contigs <- names(hdr$targets)
  true_contigs <- .get_base_name(cram_contigs)
  names(true_contigs) <- cram_contigs 

  orphans <- names(which(!true_contigs %in% rownames(spike)))
  if (!is.null(orphans)) {
    message("Orphan contigs found:")
    message(paste(orphans, collapse=", ")) 
    message("Only a subset of sequences will be represented in ", fasta)
  }

  contigs <- spike[true_contigs, "sequence"]
  names(contigs) <- names(true_contigs)
  writeXStringSet(contigs, fasta) 
  return(contigs)
    
} 


# helper fn 
.get_base_name <- function(contig_names, sep="_") {

  sapply(lapply(strsplit(contig_names, sep), `[`, 1:3), paste, collapse=sep)

}
