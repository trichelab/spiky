#' scan genomic contigs in a BAM/CRAM file
#'
#' The default workflow for spiky is roughly as follows:
#' 
#' 1. Identify and quantify the spike-in contigs in an experiment.
#' 2. Fit a model for sequence-based abundance artifacts using the spike-ins.
#' 3. Quantify raw fragment abundance on genomic contigs, and adjust per step 2.
#' 
#' scan_genomic_contigs addresses the first half of step 3. The assumption is 
#' that anything which 
#' 
#' @param bam       the BAM or CRAM file
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param param     a ScanBamParam object specifying which reads to count (NULL)
#' @param ...       additional arguments to pass to scanBamParam()
#'
#' @return          a CompressedGRangesList with bin- and spike-level coverage
#'
#' @examples
#' library(GenomicRanges)
#' data(spike, package="spiky")
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE) # swap for a CRAM 
#' res <- scan_genomic_contigs(sb, spike=spike, bins=GRanges())
#' res
#'
#' @seealso         GenomeInfoDb::keepStandardChromosomes
#' @seealso         Rsamtools::ScanBamParam
#'
#' @import          GenomicAlignments
#' @import          GenomeInfoDb
#' @import          Rsamtools
#'
#' @export
scan_genomic_contigs <- function(bam, spike, param=NULL, ...) { 

  # scan the BAM (or CRAM if supported) to determine which reads to import
  si <- seqinfo_from_header(bam)
  spike_contigs <- names(attr(find_spike_contigs(si, spike=spike), "mapping"))
  genomic_contigs <- setdiff(seqlevels(si), spike_contigs)
  if (length(genomic_contigs) == 0) {
    stop(bam, " doesn't appear to have any genomic contigs.")
  }

  # properly indexed
  bf <- BamFile(bam)

  # create appropriate filters for coverage tabulation if param=NULL
  if (is.null(param)) { 
    fl <- scanBamFlag(isDuplicate=FALSE, 
                      isPaired=TRUE, 
                      isProperPair=TRUE, ...)
    param <- ScanBamParam(flag=fl)
    bamMapqFilter(param) <- 20
  }

  # rationalize the contigs 
  gr <- as(sortSeqlevels(si[genomic_contigs]), "GRanges") # kludgey

  # restrict to only these contigs
  bamWhich(param) <- gr

  # assess coverage on these contigs (bin later)
  GenomicAlignments::coverage(BamFile(bam), param=param)

}
