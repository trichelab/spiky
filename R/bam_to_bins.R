#' create a tiled representation of a genome from the BAM/CRAM file
#' 
#' This function is meant to allow something like the following:
#' bedtools intersect -wao -a fragments.bed -b hg38_300bp_windows.bed > data.bed
#' 
#' where wao means:
#' 	-wao	Write the original A and B entries plus the number of base
#'        pairs of overlap between the two features.
#'
#' The idea is to skip the BED creation step for most runs, and just do it once.
#' In order to count reads in bins, we need bins.  
#' In order to have bins, we need to know how long the chromosomes are. 
#' In order to have a BAM or CRAM file, we need to have those same lengths.
#' This function takes advantage of all of the above to create binned ranges.
#' Note that a very recent branch of Rsamtools is required for CRAM file bins.
#' 
#' @param x       a BAM or CRAM filename (or a BamFile object)
#' @param width   the width of the bins to tile (default is 300) 
#' @param param   optional ScanBamParam (whence we attempt to extract `which`)
#' @param which   an optional GRanges restricting the bins to certain locations 
#' @param ...     additional arguments to pass on to seqinfo_from_header
#' 
#' @return        a GRangesList with y-base-pair-wide bins tiled across it
#'
#' @seealso seqinfo_from_header
#' 
#' @import IRanges
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import Rsamtools 
#' 
#' @export
bam_to_bins <- function(x, width=300, param=NULL, which=IRangesList(), ...) { 
 
  gr <- seqinfo_from_header(x, ret="gr", ...)
  if (!is.null(param)) which <- bamWhich(param) 
  if (length(which) > 0) gr <- subsetByOverlaps(gr, which)
  .renameBins(sort(unlist(tile(gr, width=width))), width=width)

}


# helper fn
.renameBins <- function(gr, width) {
  # sls <- seqlevelsStyle(gr)[1]
  # seqlevelsStyle(gr) <- "UCSC"
  names(gr) <- paste(paste0(width, "bp"), paste0("bin", seq_along(gr)), sep="_")
  # seqlevelsStyle(gr) <- sls
  return(gr)
}
