#' Binned estimation of picomoles of DNA present in cfMeDIP assays
#' 
#' This task was originally handled by the following: 
#' 
#' bedtools intersect -wao -a fragments.bed -b hg38_300bp_windows.bed > data.bed
#' 
#' The preceding translates to "Write the original A and B entries,
#' plus the number of base pairs of overlap between the two features".
#' This is essentially the same as a findOverlaps result (a Hits object),
#' but with a score consisting of the base pairs overlapped for each hit.
#' 
#' I wrote the `scan_spiked_bam` function to deal with this, and also count 
#' the number of reads that aligned to the spike-in contigs. (You're welcome.)
#' Turning the preceding per-file results into a matrix then feeds `bin_pmol`.
#' 
#' @param x     a BAM or CRAM file (the latter requires a bioc-devel fork)
#' @param bins  the bins to overlap (generated if not provided; see Details)
#' @param width width of the bins (default 300bp) to be generated if bins=NULL
#' @param mapq  minimum MAPQ value to be counted in coverage estimates (20)
#' @param ...   additional arguments (e.g. `mapq`) for computing coverage
#'
#' @return      a GRanges with average read coverage across 
#'
#' @details
#' It turns out that the fastest way to do this in Bioconductor, at least when 
#' x is a BAM or CRAM file, is to use the `bamsignals` package. For GRanges or
#' BED files, we are mostly stuck with the usual disjoin/overlap loop, but this
#' can be partly automated with bam_to_bin/bed_to_bin (using seqinfo) or just 
#' seqinfo_from_header(x). Regardless, if you want to look at a subset of bins 
#' (as opposed to all contigs or standard chromosomes), it's best to provide 
#' your own `bins` GRanges. Similarly, if you want to tally a LOT of BAM/CRAM
#' files over specific regions, it is a good idea to generate `bins` once, then
#' pass it in, although bamsignals has a convenience function that makes it
#' less important to do so (at least for speed considerations). 
#' 
#' @seealso seqinfo_from_header
#' @seealso scan_spiked_bam
#'
#' @import Rsamtools 
#' @import rtracklayer
#' @import GenomeInfoDb
#' 
#' @export
bin_pmol <- function(x, bins=NULL, width=300, mapq=20, ...) { 

  stopifnot(file.exists(x))
  res <- scan_spiked_bam(x, bins=bins, binwidth=300L, binwidth=width, ...)
  message("Computed genomic and spike coverage.")
  browser()

  # finish this modeling (call .bin_pmol_old or some such) 
  stop("This function needs to be finished") 

  # make any necessary fixes here
  return(res)

}


# helper fn
.bin_pmol_bam <- function(bam, bins=NULL, std=TRUE, ...) { 

  if (is.null(bins)) bins <- seqinfo_from_header(bam, std=std, ret="gr")
   

}  


# helper fn
.bin_pmol_gr <- function(gr, bins, std=TRUE) { 

  ABC <- disjoin(c(gr, bins))
  olA <- findOverlaps(ABC, A) 
  olB <- findOverlaps(ABC, B) 
  stop("not done yet") 

}  


### Calculates overlap with hg38 300bp windows # this is a GenomicRanges tweak


# this is just a count, can fix for submission; CRAM counts are tougher 

# big mess
# don't use until we clean this up 
bin_pmol_old <- function(...) {

  stop("this ain't done") 
  

  do.call(main, args(...))



} 

# BED file as input 
.main <- function(filename) {
  
  data <- read.table(filename, sep ="\t",  
                     header = FALSE, 
                     stringsAsFactors = FALSE)

  colnames(data) <- c("frag_chr", "frag_start", "frag_end",
                      file_path_sans_ext(basename(filename)),
                      "chr","start","end","overlap")
  
  data$frag_chr <- NULL
  data$frag_start <- NULL
  data$frag_end <- NULL
  
  ###Adjust the methylation values (overlaps/300*conc)
  window_size <- 300
  
  data$multipler <- data$overlap/window_size
  
  samples <- data[ , grepl(file_path_sans_ext(basename(filename)) , names(data))]
  adj_vals <- samples*data$multipler
  
  data <- cbind(data[,c('chr','start','end')],adj_vals)
  
  data$start <- as.factor(data$start)
  data$end <- as.factor(data$end)
  data$window <- paste0(data$chr,"_",data$start,"_", data$end)
  
  print("starting data aggregation")
  data_agg <- ddply(data,.(data$window, data$chr, data$start, data$end),numcolwise(sum), .progress = "text")
  print("finished data aggregation")
  
  print(head(data_agg))
  
  data_agg$window <- as.factor(data_agg$'data$window')
  data_agg$'data$window' <- NULL
  data_agg$chrom <- as.factor(data_agg$'data$chr')
  data_agg$'data$chr' <- NULL
  data_agg$chromStart <- as.numeric(as.character(data_agg$'data$start'))
  data_agg$'data$start' <- NULL
  data_agg$chromEnd <- as.numeric(as.character(data_agg$'data$end'))
  data_agg$'data$end' <-NULL
  print(head(data_agg))
  
  ##Now save the .bed file to compare to hg38
  data_agg$chrom <- as.character(data_agg$chrom)
  data_agg$chromStart <- as.integer(data_agg$chromStart)
  data_agg$chromEnd <- as.integer(data_agg$chromEnd)
  data_agg <- data_agg[,c("chrom", "chromStart", "chromEnd","window","adj_vals")]
  colnames(data_agg) <- c("chrom", "chromStart", "chromEnd","window",file_path_sans_ext(basename(filename)))
  print(head(data_agg))
  
  write.table(data_agg, quote = FALSE, row.names = FALSE, sep = "\t", file =  paste0("/cluster/projects/hoffmangroup/data_samanthawilson/2020_Batch3_DT/", file_path_sans_ext(basename(filename)),"_aggbywindow.bed"))
}

if (FALSE) {
  
  # this is kind of a dot-args situation -- use those
  args <- if (length(commandArgs(TRUE)) || 
              commandArgs()[length(commandArgs())] == "--args") {
    return(as.list(commandArgs(TRUE)))
  } else {
    return(list()) 
  }
} 

