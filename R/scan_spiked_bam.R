#' pretty much what it says: scan standard chroms + spike contigs from a BAM
#' 
#' @param bam       the BAM file  
#' @param mapq      minimum mapq value to count a pair (20) 
#' @param binwidth  width of the bins for chromosomal tiling (300) 
#' @param bins      a pre-tiled GRanges for binning coverage (NULL)
#' @param spike     the spike-in references, if not using the defaults (NULL) 
#' @param how       how to record spike read coverage (max or mean)? (max)
#' @param dupe      unique (FALSE), duplicte (TRUE), or all (NA) reads? (FALSE)
#' @param paired    restrict coverage to that from properly paired reads? (TRUE)
#' @param standard  restrict non-spike contigs to "standard" chromosomes? (TRUE)
#' @param ...       additional arguments to pass to scanBamFlag() 
#' 
#' @return          a CompressedGRangesList with bin- and spike-level coverage
#' 
#' @details
#'   For example (not run), one might do something like:
#'
#'   data(spike, package="spiky"); 
#'   bam <- "2021_ctl.hg38_withSpikes.bam";
#'   spiked_covg <- scan_spiked_bam(bam, mapq=20, spikes=spike);
#' 
#'   This results in a GRangesList object with 300bp-binned coverage on the 
#'   standard (chr1-22, chrX, chrY, chrM) chromosomes (as determined by the
#'   GenomeInfoDb::standardChromosomes() function against the assembly defined
#'   in the BAM or CRAM file, by default; if desired, a user can scan all 
#'   genomic contigs by setting standard=FALSE when calling the function). 
#'   By default, the mean base-level coverage of genomic bins is reported, 
#'   and the maximum spike-level coverage is reported, though this can also 
#'   be adjusted as needed. The results then inform the reliability of
#'   measurements from replicate samples in multiple labs, as well as the 
#'   adjusted quantitative coverage in each bin once the absolute quantity
#'   of captured cell-free methylated DNA has been fit by model_glm_pmol and 
#'   predict_pmol. In some sense, this function converts BAMs/CRAMs into usable 
#'   data structures for high-throughput standardized cfMeDIP experiments.
#' 
#' @seealso         GenomeInfoDb::keepStandardChromosomes
#' @seealso         Rsamtools::ScanBamParam
#' 
#' @import          GenomicAlignments
#' @import          GenomeInfoDb
#' @import          Rsamtools
#' 
#' @export
scan_spiked_bam <- function(bam, mapq=20, binwidth=300L, bins=NULL, spike=NULL, how=c("max", "mean"), dupe=FALSE, paired=TRUE, standard=TRUE, ...) { 

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

  # properly indexed
  bf <- BamFile(bam)

  # get a GRanges corresponding to our contig selections for the BAM/CRAM file
  gr <- get_merged_gr(si=si, standard=standard) # scanBam positive selections

  # tabulate spike coverage onto this subset of the filter GRanges from above
  spike_gr <- subset(gr, genome(gr)[as(seqnames(gr), "character")] == "spike")

  # create appropriate filters for coverage tabulation
  fl <- scanBamFlag(isDuplicate=dupe, isPaired=paired, isProperPair=paired, ...)
  bp <- ScanBamParam(flag=fl, which=gr)
  bamMapqFilter(bp) <- mapq

  # assess coverage on contigs we care about (bin the Rles later)
  covg <- get_spiked_coverage(bf=bf, bp=bp, gr=gr)

  # bins for coverage
  if (is.null(bins)) bins <- tile_bins(gr=gr, binwidth=binwidth)

  # genomic coverage is averaged across each bin
  binned <- get_binned_coverage(bins, covg)

  # spikes get summarized as max or mean across each
  spiked <- get_spike_depth(covg, spike_gr, how=how, spike=spike) 
    
  # toss them into a list of summary results
  res <- GRangesList(genomic=binned, spikes=spiked) # redundant mcols :-( 
  return(res)

} 
