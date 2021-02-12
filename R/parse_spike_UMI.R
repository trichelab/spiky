#' parse out the forward and reverse UMIs and contig for a BED/BAM
#'
#' @param   UMI   a vector of UMIs
#' @param   pos   optional vector of positions (else all are set to 1) 
#' @param   seqs  optional vector of read sequences (else widths default to 96)
#' 
#' @return  a GRanges
#' 
#' @import Biostrings
#' @import GenomicRanges
#' 
#' @export 
parse_spike_UMI <- function(UMI, pos=NULL, seqs=NULL) { 
  nch <- nchar(UMI)
  UMIs <- strsplit(substr(UMI, nch - 10, nch), "_")
  x <- data.frame(UMI1 = DNAStringSet(vapply(UMIs, `[`, character(1), 1)),
                  UMI2 = DNAStringSet(vapply(UMIs, `[`, character(1), 2)),
                  chrom = substr(UMI, 1, nch - 11))

  x$chromStart <- 1
  if (!is.null(pos)) x$chromStart <- pos
  widths <- rep(96, nrow(x))
  if (!is.null(seqs)) x$seq <- DNAStringSet(seqs)
  if (!is.null(seqs)) widths <- width(x$seq)
  x$chromEnd <- x$chromStart + widths - 1 
  x$name <- UMI 
  x$score <- 1 # subject 6547 

  return(as(x, "GRanges"))
} 
