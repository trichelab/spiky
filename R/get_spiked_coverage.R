#' tabulate coverage across assembly and spike contig subset in natural order
#'
#' Refactored from scan_spiked_bam, this is a very simple wrapper
#'
#' @param bf  the BamFile object 
#' @param bp  the ScanBamParam object 
#' @param gr  the GRanges with sorted seqlevels
#' 
#' @return    a list of Rles 
#' 
#' @seealso   scan_spiked_bam
#' @seealso   coverage
#' 
#' @examples
#' sb <- system.file("extdata", "example.spike.bam", package="spiky", 
#'                   mustWork=TRUE)
#' si <- seqinfo_from_header(sb) 
#' genome(si) <- "spike"
#' mgr <- get_merged_gr(si)
#'
#' fl <- scanBamFlag(isDuplicate=FALSE, isPaired=TRUE, isProperPair=TRUE)
#' bp <- ScanBamParam(flag=fl)
#' bamMapqFilter(bp) <- 20
#' 
#' get_spiked_coverage(sb, bp=bp, gr=mgr)
#'
#' @import    GenomicAlignments
#' 
#' @export
get_spiked_coverage <- function(bf, bp, gr) { 
  
  message("Tabulating read pair coverage (may take a while)...", appendLF=FALSE)
  covg <- coverage(bf, param=bp)[seqlevels(gr)]
  message("Done.") 
  return(covg)

}
