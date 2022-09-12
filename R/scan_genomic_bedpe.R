#' Scan genomic BEDPE
#'
#'
#' @param bedpe       the BEDPE file path, or output from read_bedpe()
#' @param bin       Bin reads? (TRUE)
#' @param binwidth  width of the bins for chromosomal tiling (300)
#' @param bins      a pre-tiled GRanges for binning coverage (NULL)
#' @param standard  restrict non-spike contigs to "standard" chromosomes? (TRUE)
#' @param genome    Name of genome (default hg38)
#'
#' @return          a GRanges with coverage
#'
#' @examples
#'
#' fl <- system.file("extdata", "example_chr21_bedpe.bed.gz", package="spiky",mustWork=TRUE)
#' scan_genomic_bedpe(fl) # will warn user about spike contigs
#'
#'
#' @export
scan_genomic_bedpe <- function(bedpe, bin=TRUE, binwidth=300L, bins=NULL, standard=TRUE,genome="hg38"){
  if (!is(bedpe, "Pairs")) {bedpe <- read_bedpe(bedpe)}
  genomic <- convertPairedGRtoGR(bedpe)

  g_covg <- coverage(genomic)

  g_si <- seqinfo(g_covg)
  seqlengths(g_si) <- sapply(g_covg, length)
  if (standard) seqlevels(g_si) <- seqlevels(keepStandardChromosomes(g_si))
  genome(g_si) <- genome
  g_gr <- as(g_si,"GRanges")
  if (is.null(bins)) bins <- tile_bins(gr=g_gr, binwidth=binwidth)

  if (bin){
    # genomic coverage is averaged across each bin
    if (length(bins) > 0) {
      message("Binning genomic coverage...")
      binned <- get_binned_coverage(bins, g_covg)
      message("Done.")
    } else {
      message("Empty bins provided, skipping genomic binning.")
      binned <- bins
    }
  } else {binned <- bins}

  gr <- as(binned,"GRanges")
  names(gr) <- seqnames(gr)
  return(gr)

}


