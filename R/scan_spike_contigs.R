#' pretty much what it says: scan spike contigs from a BAM or CRAM file
#'
#' default workflow is
#'
#' 1. scan spike contigs and count fragments per contig or per bin.
#' 2. fit the appropriate model for adjusting genomic contigs based on spikes.
#' 3. scan and adjust binned fragment tallies along genomic contigs per above.
#'
#' scan_spike_contigs implements step 1.

#'
#' @param bam       the BAM or CRAM filename, or a vector of such filenames
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param how       how to summarize the per-spike coverage (max)
#' @param param     a ScanBamParam object, or NULL (will default to MAPQ=20 etc)
#' @param ...       additional arguments to pass to scanBamFlag()
#' @param mc.cores  Number of cores to run on (default 16)
#'
#' @details
#' If multiple BAM or CRAM filenames are provided, all indices will be
#' checked before attempting to run through any of the files.
#'
#' @return          a CompressedGRangesList with bin- and spike-level coverage
#'
#' @examples
#' library(GenomicRanges)
#' data(spike, package="spiky")
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE) # switch to a CRAM
#' res <- scan_spike_contigs(sb, spike=spike) # use default ScanBamParam
#' summary(res)
#'
#' @seealso         Rsamtools::ScanBamParam
#'
#' @import          GenomicAlignments
#' @import          GenomeInfoDb
#' @import          Rsamtools
#'
#' @export
scan_spike_contigs <- function(bam, spike, how="max", param=NULL, mc.cores=16,...) {

  # can be smoother but:
  if (length(bam) > 1) {
    indices <- sub("bam$", "bam.bai", bam)
    indices <- sub("cram$", "cram.crai", bam)
    if (!all(file.exists(indices))) {
      missed <- indices[!file.exists(indices)]
      stop("Missing index files: ", paste(missed, collapse=", "))
    } else {
      if (is.null(names(bam))) names(bam) <- bam
      return(lapply(bam, scan_spike_contigs, spike=spike))
    }
  }

  # scan the BAM (or CRAM if supported) to determine which reads to import
  si <- seqinfo_from_header(bam)
  spikes <- find_spike_contigs(si, spike=spike)
  if (length(spikes) > 0) {
    mappings <- attr(spikes, "mappings")
    attr(spikes, "mappings") <- NULL
    genome(si)[spikes] <- "spike"
  }

  # create appropriate filters for coverage tabulation if param=NULL
  if (is.null(param)) {
    fl <- scanBamFlag(isDuplicate=FALSE,
                      isPaired=TRUE, ...)
    param <- ScanBamParam(flag=fl)
    bamMapqFilter(param) <- 20
  }

  orig_spike_contigs <- subset(seqlevels(si), genome(si) == "spike")
  if (length(orig_spike_contigs) == 0) {
    # empty coverage list
    s <- as(S4Vectors::SimpleList(), "SimpleRleList")
  } else {
    # rationalize the contigs but don't override user supplied bamWhich
    gr <- as(sortSeqlevels(si[orig_spike_contigs]), "GRanges") # kludgey
    if (length(bamWhich(param)) == 0) bamWhich(param) <- gr
    s <- GenomicAlignments::coverage(BamFile(bam), param=param)
  }
  s_depth <- get_spike_depth(s,spike=spike,how=how)
  return(as(s_depth,"GRanges"))
}
