#' for BAM/CRAM files with renamed contigs, we need to rename `spike` rows
#' 
#' This function does that. 
#' 
#' @param x         a BAM/CRAM file, hopefully with an index
#' @param spike     a DataFrame where spike$sequence is a DNAStringSet
#' 
#' @return          a DataFrame with renamed contigs (rows)
#' 
#' @seealso generate_spike_fasta
#' 
#' @import Rsamtools
#' @import Biostrings
#' 
#' @export
rename_spikes <- function(x, spike) { 

  cram <- BamFile(x) 
  hdr <- scanBamHeader(cram) 
  cram_contigs <- names(hdr$targets)
  true_contigs <- .get_base_name(cram_contigs)
  names(true_contigs) <- cram_contigs 
  
  orphans <- names(which(!true_contigs %in% rownames(spike)))
  if (!is.null(orphans)) {
    message("Orphan contigs found:")
    message(paste(orphans, collapse=", ")) 
    message("Only a subset of sequences will be represented in ", fasta)
  }

  newspikes <- spike[true_contigs, ] 
  rownames(newspikes) <- cram_contigs 
  names(newspikes$sequence) <- cram_contigs
  return(newspikes) 

} 


# helper fn 
.get_base_name <- function(contig_names, sep="_") {

  sapply(lapply(strsplit(contig_names, sep), `[`, 1:3), paste, collapse=sep)

}
