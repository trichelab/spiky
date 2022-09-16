#' Convert Pairs to GRanges
#'
#' @param pairs          the Pairs object
#'
#' @return            a GRanges
#'
#' @export
convertPairedGRtoGR <- function(pairs){
  ir <- IRanges(start=start(first(pairs)),end=end(second(pairs)))
  gr <- GRanges(seqnames=seqnames(first(pairs)),ranges=ir)
  return(gr)
}
