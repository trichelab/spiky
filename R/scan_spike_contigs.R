#' pretty much what it says: scan standard chroms + spike contigs from a BAM
#'
#' Note: behind the scenes, this is being refactored into scan_spike_contigs 
#' and scan_genomic_contigs. Once that is done, perhaps before release, the 
#' default workflow will switch to 
#' 
#' 1. scan spike contigs and count fragments per contig or per bin.
#' 2. fit the appropriate model for adjusting genomic contigs based on spikes.
#' 3. scan and adjust binned fragment tallies along genomic contigs per above.
#' 
#' This approach decouples binning schemes from model generation (using spikes) 
#' and model-based adjustment (using genomic fragment counts), decreasing code
#' complexity while increasing the opportunities for caching & parallelization.
#' 
#' @param bam       the BAM file
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param param     a ScanBamParam object, or NULL (will default to MAPQ=20 etc)
#' @param ...       additional arguments to pass to scanBamFlag()
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
#' @details
#'   add CRAM example here -- tested & works with reheadered spike CRAMs.
#'   Slower than one might like however.
#' 
#' @seealso         Rsamtools::ScanBamParam
#'
#' @import          GenomicAlignments
#' @import          GenomeInfoDb
#' @import          Rsamtools
#'
#' @export 
scan_spike_contigs <- function(bam, spike, param=NULL, ...) {

  # scan the BAM (or CRAM if supported) to determine which reads to import
  si <- seqinfo_from_header(bam)
  spikes <- find_spike_contigs(si, spike=spike)
  if (length(spikes) > 0) {
    mappings <- attr(spikes, "mappings")
    attr(spikes, "mappings") <- NULL
    genome(si)[spikes] <- "spike"
  } else {
    stop(bam, " doesn't appear to have any spike-ins among its contigs.")
  }

  # create appropriate filters for coverage tabulation if param=NULL
  if (is.null(param)) { 
    fl <- scanBamFlag(isDuplicate=FALSE, 
                      isPaired=TRUE, 
                      isProperPair=TRUE, ...)
    param <- ScanBamParam(flag=fl)
    bamMapqFilter(param) <- 20
  }

  # rationalize the contigs 
  orig_spike_contigs <- subset(seqlevels(si), genome(si) == "spike")
  gr <- as(sortSeqlevels(si[orig_spike_contigs]), "GRanges") # kludgey

  # restrict to only the spike contigs
  bamWhich(param) <- gr

  # assess coverage on spike contigs (bin later)
  GenomicAlignments::coverage(BamFile(bam), param=param)

}
