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
#' Using other versions of Rsamtools will yield an error on CRAM files.
#' 
#' @param x         a CRAM file, hopefully with an index
#' @param fasta     the filename for the resulting FASTA ("spikes.fa")
#' 
#' @return          invisibly, a DNAStringSet as exported to `fasta` 
#'
#' @seealso         rename_contigs
#' 
#' @import          Biostrings 
#' @import          Rsamtools 
#' 
#' @export
generate_spike_fasta <- function(x, fasta="spike_contigs.fa") { 

  cram <- x 
  if (!is(cram, "BamFile")) cram <- BamFile(cram) 
  hdr <- scanBamHeader(cram) 
  cram_contigs <- names(hdr$targets)
  
  data(spike, package="spiky") 
  newspikes <- rename_spikes(cram, spike=spike)

  contigs <- newspikes[cram_contigs, "sequence"]
  names(contigs) <- cram_contigs
  writeXStringSet(contigs, fasta)
  message("Wrote ", fasta) 
  invisible(contigs)
   
} 


# helper fn 
.get_base_name <- function(contig_names, sep="_") {

  sapply(lapply(strsplit(contig_names, sep), `[`, 1:3), paste, collapse=sep)

}
