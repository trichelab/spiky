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
#' @param bam       the BAM or CRAM filename, or a vector of them
#' @param spike     the spike-in reference database (e.g. data(spike))
#' @param param     a ScanBamParam object specifying which reads to count (NULL)
#' @param ...       additional arguments to pass to scanBamFlag()
#'
#' @return          a CompressedGRangesList with bin- and spike-level coverage
#'
#' @details
#' If multiple BAM or CRAM filenames are provided, all indices will be 
#' checked before attempting to run through any of the files. 
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
scan_genomic_contigs <- function(bam, spike, param=NULL, ...) {

  # can be smoother but:
  if (length(bam) > 1) {
    indices <- sub("bam$", "bam.bai", bam)
    indices <- sub("cram$", "cram.crai", bam)
    if (!all(file.exists(indices))) {
      missed <- indices[!file.exists(indices)]
      stop("Missing index files: ", paste(missed, collapse=", "))
    } else {
      if (is.null(names(bam))) names(bam) <- bam
      return(lapply(bam, scan_genomic_contigs, spike=spike))
    }
  }

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

  # assess coverage on these contigs (bin later)
  return(GenomicAlignments::coverage(BamFile(bam), param=param))

}
