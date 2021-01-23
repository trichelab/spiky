#' get a GRanges of (by default, standard) chromosomes from seqinfo
#' 
#' refactored from scan_spiked_bam to clarify information flow
#' 
#' @param si        seqinfo, usually from a BAM/CRAM file with spike contigs
#' @param standard  trim to standard chromosomes? (TRUE) 
#' 
#' @return          GRanges with two genomes: the organism assembly and "spike"
#'
#' @import          GenomeInfoDb
#' 
#' @export
get_merged_gr <- function(si, standard=TRUE) { 

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
  names(gr) <- as.character(seqnames(rename_spike_seqlevels(gr, spike=spike)))
  return(gr)

} 
