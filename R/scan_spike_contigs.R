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
#' @param bam_files       the BAM or CRAM file, or list of BAMs/CRAMs with the same header
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param param     a ScanBamParam object, or NULL (will default to MAPQ=20 etc)
#' @param mc.cores  the number of cores to use (default is minimum of 16 and number of BAMs)
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
scan_spike_contigs <- function(bam_files, spike, param=NULL, mc.cores=16,...) {

  # Grab the first bam in the list
  if (is.list(bam_files)){bam <- unlist(bam_files[1])}else{ bam<-bam_files }

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
                      isPaired=TRUE,
                      isProperPair=TRUE, ...)
    param <- ScanBamParam(flag=fl)
    bamMapqFilter(param) <- 20
  }

  orig_spike_contigs <- subset(seqlevels(si), genome(si) == "spike")
  if (length(orig_spike_contigs) == 0) {
    # empty coverage list
    return(as(S4Vectors::SimpleList(), "SimpleRleList"))
  } else {
    # rationalize the contigs but don't override user supplied bamWhich
    gr <- as(sortSeqlevels(si[orig_spike_contigs]), "GRanges") # kludgey
    if (length(bamWhich(param)) == 0) bamWhich(param) <- gr

    # number of cores to use
    mc.cores <- min(mc.cores, length(bam_files))
    # Multi-thread apply computes coverage on all bam files.
    return(mclapply(bam_files,FUN=function(x)
      {GenomicAlignments::coverage(BamFile(x), param=param)},mc.cores=mc.cores))
  }
}
