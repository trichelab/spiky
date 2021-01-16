#' pretty much what it says: scan standard chroms + spike contigs from a BAM
#' 
#' @param bam       the BAM file  
#' @param mapq      minimum mapq value to count a pair (20) 
#' @param binwidth  width of the bins for chromosomal tiling (300) 
#' @param bins      a pre-tiled GRanges for binning coverage (NULL)
#' @param spike     the spike-in references, if not using the defaults (NULL) 
#' @param ...       additional arguments to pass to ScanBamParam() 
#' 
#' @return          base-level read pair coverage on chroms + spikes in the BAM
#' 
#' @details
#'   For example (not run), one might do something like:
#'
#'   data(spike, package="spiky")
#'   bam <- "2021_ctl.hg38_withSpikes.bam" # a merged sorted BAM
#'   spiked_covg <- scan_spiked_bam(bam, mapq=20, spikes=spike)
#' 
#' @seealso         GenomeInfoDb::keepStandardChromosomes
#' @seealso         Rsamtools::ScanBamParam
#' 
#' @import          GenomicAlignments
#' @import          GenomeInfoDb
#' @import          Rsamtools
#' 
#' @export
scan_spiked_bam <- function(bam, mapq=20, binwidth=300L, bins=NULL, spike=NULL, ...) { 

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

  # currently we drop nonstandard contigs for this process
  standard_chroms <- seqlevels(keepStandardChromosomes(si))
  orig_spike_contigs <- subset(seqlevels(si), genome(si) == "spike")
  merged_contigs <- sortSeqlevels(si[c(standard_chroms, orig_spike_contigs)])
  gr <- as(merged_contigs, "GRanges")

  # use standardized element names, but keep the seqlevels found in the BAM:
  names(gr) <- as.character(seqnames(rename_spike_seqlevels(gr, spike=spike)))
  bp <- ScanBamParam(which=gr, ...)
  bamMapqFilter(bp) <- mapq
  bf <- BamFile(bam)

  # the fact I have to include like this suggests optimization is called for
  message("Tabulating read pair coverage (may take a while)...", appendLF=FALSE)
  system.time(cvg <- coverage(bf, param=bp)[seqlevels(gr)])
  message("Done.") 

  # welp, here's your unbinned coverage -- may be easiest to just bin the Rles?
  if (is.null(bins)) { 
    std <- seqinfo(gr)[standard_chroms]
    message("Tiling ", binwidth, "bp bins across the genome...", appendLF=FALSE)
    system.time(bins <- tileGenome(std, tilewidth=binwidth, 
                                   cut.last.tile.in.chrom=TRUE))
    message("Done.") 
  }

  # genomic coverage is averaged across each bin
  message("Binning genomic coverage...", appendLF=FALSE)
  system.time(binned <- binnedAverage(bins, numvar=cvg[seqlevels(bins)], 
                                      varname="coverage"))
  binned <- split(binned, seqnames(binned))
  message("Done.") 
 
  # spikes get summarized as median coverage across each spike
  message("Summarizing spike-in counts...", appendLF=FALSE)
  system.time(spiked <- sapply(cvg[orig_spike_contigs], median))
  names(spiked) <- get_base_name(names(spiked))
  message("Done.") 

  # toss them into a list along with a QC of max bin covg per chrom 
  maxcovg <- sapply(binned, function(x) max(x$coverage))
  res <- list(genomic=binned, spikes=spiked, maxcovg=maxcovg)
  return(res)

} 
