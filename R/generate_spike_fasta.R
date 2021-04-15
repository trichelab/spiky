#' for CRAM files, a FASTA reference is required to decode; this builds that
#'
#' A FASTA reference is *not* always needed, so long as .crai indices
#' are available for all contigs in the CRAM. See `spike_counts` for a fast
#' and convenient alternative that extracts spike coverage from index stats.
#' However, spike_counts has its own issues, and it's better to use fragments.
#'
#' If the contigs in a CRAM have even slightly different names from those in
#' the reference, decoding will fail.  In some cases there are multiple names
#' for a given contig (which raises the question of whether to condense them),
#' and thus the same reference sequence decodes multiple contig names.
#'
#' This function generates an appropriate spike reference for a BAM or CRAM,
#' using BAM/CRAM headers to figure out which references are used for which.
#'
#' At the moment, CRAM support in Rsamtools only exists in the GitHub branch:
#'
#' BiocManager::install("Bioconductor/Rsamtools@cram")
#'
#' Using other versions of Rsamtools will yield an error on CRAM files.
#'
#' Note that for merged genomic + spike reference BAMs/CRAMs, this function
#' will only attempt to generate a FASTA for the spike contigs, not reference.
#' If your reference contigs are screwed up, talk to your sequencing people,
#' and keep better track of the FASTA reference against which you compress!
#'
#' @param bam       a BAM or CRAM file, hopefully with an index
#' @param spike     the spike contig database (mandatory as of 0.9.99)
#' @param assembly  optional BSgenome or seqinfo with reference contigs (NULL)
#' @param fa        the filename for the resulting FASTA ("spikes.fa")
#'
#' @return          invisibly, a DNAStringSet as exported to `fa`
#'
#' @examples
#'
#' library(GenomicRanges)
#' data(spike, package="spiky")
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' outFasta <- paste(system.file("extdata", package="spiky", mustWork=TRUE),"/spike_contigs.fa",sep="")
#' show(generate_spike_fasta(sb, spike=spike,fa=outFasta))
#'
#' @seealso         rename_contigs
#'
#' @import          Biostrings
#' @import          Rsamtools
#'
#' @export
generate_spike_fasta <- function(bam, spike, assembly=NULL, fa="spike_contigs.fa"){

  if (!is(bam, "BamFile")) bam <- BamFile(bam)
  hdr <- scanBamHeader(bam)
  bam_contigs <- names(hdr$targets)
  if (!is.null(assembly)) {
    bam_contigs <- setdiff(bam_contigs, seqlevels(assembly))
  }

  # are there spikes in it?
  if (length(bam_contigs) == 0) stop("No spike contigs found in ", bam, "!")
  newspikes <- rename_spikes(bam, spike=spike) # this is probably dumb

  contigs <- newspikes[bam_contigs, "sequence"]
  names(contigs) <- bam_contigs
  writeXStringSet(contigs, fa)
  message("Wrote ", fa)
  invisible(contigs)

}


# helper fn
.get_base_name <- function(contig_names, sep="_") {

  vapply(lapply(strsplit(contig_names, sep), `[`, seq_len(3)),
         paste, character(1), collapse=sep)

}
