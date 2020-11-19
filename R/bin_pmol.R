#' This function replaces the following line in bedtools:
#' 
#' bedtools intersect -wao -a fragments.bed -b hg38_300bp_windows.bed > data.bed
#' 
#' Which translates to "Write the original A and B entries,
#' plus the number of base pairs of overlap between the two features".
#' 
#' Which translates to 
### Calculates overlap with hg38 300bp windows # this is a GenomicRanges tweak


# this is just a count, can fix for submission; CRAM counts are tougher 

# big mess
# don't use until we clean this up 
bin_pmol <- function(...) {

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

