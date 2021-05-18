#' scan genomic contigs in a BAM/CRAM file
#'
#' The default workflow for spiky is roughly as follows:
#'
#' 1. Identify and quantify the spike-in contigs in an experiment.
#' 2. Fit a model for sequence-based abundance artifacts using the spike-ins.
#' 3. Quantify raw fragment abundance on genomic contigs, and adjust per step 2.
#'
#' scan_genomic_contigs addresses the first half of step 3. The assumption is
#' that anything which isn't a spike contig, is a genomic contig.  This isn't
#' necessarily true, so the user can also supply a ScanBamParam object for the
#' `param` argument and restrict scanning to whatever contigs they wish, which
#' also allows for non-default MAPQ, pairing, and quality filters.
#'
#' @param bam_files       the BAM or CRAM file, or list of BAMs/CRAMs with the same header
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param param     a ScanBamParam object specifying which reads to count (NULL)
#' @param mc.cores  the number of cores to use (default is minimum of 16 and number of BAMs)
#' @param ...       additional arguments to pass to scanBamFlag()
#'
#' @return          a CompressedGRangesList with bin- and spike-level coverage
#'
#' @examples
#'
#' library(Rsamtools)
#' data(spike, package="spiky")
#'
#' fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
#'                   mustWork=TRUE)
#' scan_genomic_contigs(fl, spike=spike) # will warn user about spike contigs
#'
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' scan_genomic_contigs(sb, spike=spike) # will warn user about genomic contigs
#'
#' @seealso         Rsamtools::ScanBamParam
#' @import          GenomicAlignments
#' @import          Rsamtools
#'
#' @export
scan_genomic_contigs <- function(bam_files, spike, param=NULL, ...) {

  # Grab the first bam in the list
  if (is.list(bam_files)){bam <- unlist(bam_files[1])}else{ bam<-bam_files }

  # scan the BAM (or CRAM if supported) to determine which reads to import
  si <- seqinfo_from_header(bam)
  mappings <- attr(find_spike_contigs(si, spike=spike), "mapping")
  spike_contigs <- names(mappings)
  genomic_contigs <- seqlevels(si)
  if (length(spike_contigs) > 0) {
    genomic_contigs <- setdiff(genomic_contigs, spike_contigs)
  }

  if (length(genomic_contigs) == 0) {
    # empty coverage list
    warning(bam, " doesn't appear to have any genomic contigs.")
    return(as(S4Vectors::SimpleList(), "SimpleRleList"))
  }

  bf <- BamFile(bam)
  if (is.null(param)) {
    fl <- scanBamFlag(isDuplicate=FALSE,
                      isPaired=TRUE,
                      isProperPair=TRUE, ...)
    param <- ScanBamParam(flag=fl)
    bamMapqFilter(param) <- 20
  }

  # rationalize the contigs (but do not replace user-supplied Ranges)
  gr <- as(sortSeqlevels(si[genomic_contigs]), "GRanges") # kludgey
  if (length(bamWhich(param)) == 0) bamWhich(param) <- gr

  # number of cores to use
  mc.cores <- min(mc.cores, length(bam_files))

  # assess coverage on these contigs (bin later)
  return(mclapply(bam_files,FUN=function(x)
  {GenomicAlignments::coverage(BamFile(x), param=param)},mc.cores=mc.cores))

}
