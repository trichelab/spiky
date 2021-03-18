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
#' BiocManager::install("Bioconductor/Rhtslib@enable-cram")
#' BiocManager::install("Bioconductor/Rsamtools@cram") 
#'
#' Using other versions of Rsamtools will yield an error on CRAM files.
#' 
#' Note that for merged genomic + spike reference BAMs/CRAMs, this function 
#' will only attempt to generate a FASTA for the spike contigs, not reference.
#' If your reference contigs are screwed up, talk to your sequencing people, 
#' and keep better track of the FASTA reference against which you compress!
#' 
#' @param x         a CRAM file, hopefully with an index
#' @param assembly  optional BSgenome or seqinfo with reference contigs (NULL)
#' @param fasta     the filename for the resulting FASTA ("spikes.fa")
#' 
#' @return          invisibly, a DNAStringSet as exported to `fasta` 
#'
#' @examples 
#'
#'  
#' @seealso         rename_contigs
#' 
#' @import          Biostrings 
#' @import          Rsamtools 
#' 
#' @export
generate_spike_fasta <- function(x, assembly=NULL , fasta="spike_contigs.fa") { 

  cram <- x 
  if (!is(cram, "BamFile")) cram <- BamFile(cram) 
  hdr <- scanBamHeader(cram) 
  cram_contigs <- names(hdr$targets)
 
  if (!is.null(assembly)) {
    cram_contigs <- setdiff(cram_contigs, seqlevels(assembly))
  }
  
  # are there spikes in it? 
  if (length(cram_contigs) == 0) stop("No spike contigs found in ", x, "!")

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

  vapply(lapply(strsplit(contig_names, sep), `[`, seq_len(3)), 
         paste, character(1), collapse=sep)

}
