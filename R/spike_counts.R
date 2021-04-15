#' use the index of a spiked BAM/CRAM file for spike contig coverage 
#' 
#' It dawned on me one day that we don't even have to bother reading the file
#' if we have an index for a spiked BAM/CRAM result, since any fragments that 
#' map properly to the spike contigs are generated from synthetic templates.
#' This function takes an index and a spike database (usually a DataFrame) as
#' inputs and provides a rough coverage estimate over "rehabilitated" contig
#' names (i.e., canonicalized contigs mapping to the database) as its output. 
#' 
#' @param     bam             the BAM or CRAM file (MUST HAVE AN INDEX) 
#' @param     spike           a data.frame, DataFrame, or similar with spikes
#' @param     sep             separator character in contig names ("_")
#' @param     ref             reference name for spike genome ("spike")
#' @param     verbose         be verbose? (FALSE) 
#' @param     dump_idx        dump the renamed idxstats to aggregate? (FALSE)
#' 
#' @examples
#' data(spike, package="spiky")
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE) 
#' spike_counts(sb, spike=spike)
#' 
#' @details
#' The argument `spike` has no default since we are attempting to refactor the
#' spike-in databases into their own data packages and allow more general use.
#' 
#' @return                    a GRanges of spike-in contig read counts
#'
#' @import GenomicRanges
#' @import Rsamtools
#' 
#' @export
spike_counts <- function(bam, spike, sep="_", ref="spike", verbose=FALSE, dump_idx=FALSE) {

  if (!is(bam, "BamFile")) bam <- BamFile(bam)
  if (!file.exists(bam$index)) {
    message("You need a valid index for ", bam$path, " to use this function.")
    message("(Please do NOT use Rsamtools::indexBam() if it is a CRAM file!)")
    stop(bam$index, " was not found; aborting.")
  }
  if (!all(c("sequence", "methylated", "CpGs") %in% names(spike)) |
      !all(grepl(sep, rownames(spike)))) {
    stop("Your spike-in database is missing either columns or names; aborting.")
  }

  if (verbose) message("Scanning ", bam$index, "...")
  isb <- idxstatsBam(bam) 
  sbn <- get_base_name(isb[, "seqnames"])
  spike_contigs <- sbn[which(sbn %in% rownames(spike))]
  isb <- subset(isb, seqnames %in% names(spike_contigs))
  isb$seqnames <- names(spike_contigs)
  rownames(isb) <- spike_contigs
  isb$methylated <- spike[rownames(isb), "methylated"] 
  isb$id <- gsub("\\.(bam|cram)", "", basename(bam$path))
  if (verbose) message(length(spike_contigs), " spike-in contigs identified.")

  # used by scan_spike_counts
  if (dump_idx) {
    return(isb)
  } else { 
    # prep seqinfo
    isbsi <- as(isb, "Seqinfo") 
    seqlengths(isbsi) <- isb$seqlength
    genome(isbsi) <- ref

    # prep covg
    isb$start <- 1
    isb$end <- isb$seqlength
    isb$unmapped <- NULL # irrelevant
    names(isb) <- sub("mapped", "coverage", names(isb))
    covg <- as(isb[, c("seqnames", "start", "end", "coverage", "methylated")], 
               "GRanges")
    seqinfo(covg) <- isbsi[seqlevels(covg)]

    # done
    return(covg)
  }

}
