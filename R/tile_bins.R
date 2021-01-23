#' Tile the assembly-based contigs of a merged assembly/spike GRanges.
#'
#' refactored out of scan_spiked_bam for more explicit information flow
#' 
#' @param gr          the GRanges
#' @param binwidth    bin width to tile (default is 300)
#' 
#' @return            a GRanges of bins
#'
#' @export
tile_bins <- function(gr, binwidth=300L) {
 
  chroms <- subset(seqlevels(gr), genome(gr)[seqlevels(gr)] != "spike")
  si <- seqinfo(gr)[chroms] 
  message("Tiling ", binwidth, "bp bins across the genome...", appendLF=FALSE)
  bins <- tileGenome(si, tilewidth=binwidth, cut.last.tile.in.chrom=TRUE)
  message("Done.") 
  return(bins)

}
