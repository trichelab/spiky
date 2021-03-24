#' get a GRanges of (by default, standard) chromosomes from seqinfo
#' 
#' refactored from scan_spiked_bam to clarify information flow
#' 
#' @param si        seqinfo, usually from a BAM/CRAM file with spike contigs
#' @param standard  trim to standard chromosomes? (TRUE) 
#' @param spike     database of spike-in standard sequence features (spike)
#' 
#' @return          GRanges with two genomes: the organism assembly and "spike"
#'
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky", 
#'                   mustWork=TRUE) 
#' si <- seqinfo_from_header(sb) 
#' genome(si) <- "spike" # no genomic contigs
#' data(spike) # will soon be required anyhow 
#' get_merged_gr(si, spike=spike) # note canonicalized spikes
#' 
#' @details
#' By default, get_merged_gr will return a GRanges with "standardized" 
#' genomic and spike contig names (i.e. genomic chr1-22, X, Y, M, and 
#' the canonical spike names in data(spike, package="spiky")).
#' 
#' The constraint to "standard" chromosomes on genomic contigs can be 
#' removed by setting `standard` to FALSE in the function arguments.
#'
#' @import          GenomeInfoDb
#' 
#' @export
get_merged_gr <- function(si, standard=TRUE, spike=NULL) { 

  # assembly contigs 
  if (standard) {
    chrom_contigs <- seqlevels(keepStandardChromosomes(si))
  } else {
    chrom_contigs <- seqlevels(si)
  }

  # spike contigs using the names found in the merged BAM/CRAM file
  orig_spike_contigs <- subset(seqlevels(si), genome(si) == "spike")

  # the collection of contigs we want to keep for downstream analysis
  merged_contigs <- sortSeqlevels(si[c(chrom_contigs, orig_spike_contigs)])

  # coerce this to a GRanges object 
  gr <- as(merged_contigs, "GRanges")

  # use standardized element names, but keep the seqlevels found in the BAM:
  if (is.null(spike)) data(spike)
  names(gr) <- as.character(seqnames(rename_spike_seqlevels(gr, spike=spike)))
  return(gr)

} 
