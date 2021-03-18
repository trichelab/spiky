#' run spike_counts on BAM/CRAM files and shape the results for model_glm_pmol
#' 
#' Typically one will want to fit a correction model to multiple samples. 
#' This function eases this task by merging the output of spike_counts into 
#' a data.frame that model_glm_pmol can directly fit. 
#' 
#' @param files   a vector of BAM/CRAM file names
#' @param spike   a spike-in database 
#' @param sep     the separator for spike-in contig names ("_") 
#' 
#' @return        a data.frame with columns "frag_grp", "id", and "read_count"
#'
#' @examples
#' data(spike) 
#' library(GenomicRanges)
#' sb <- system.file("extdata", "example.spike.bam", package="spiky", 
#'                   mustWork=TRUE)
#' scan_spike_counts(sb, spike=spike)
#' fit <- model_glm_pmol(scan_spike_counts(sb, spike=spike)) 
#'
#' @export
scan_spike_counts <- function(files, spike, sep="_") { 

  stopifnot(all(file.exists(files)))
  ssc <- lapply(files, spike_counts, spike=spike, dump_idx=TRUE)
  do.call(rbind, lapply(lapply(ssc, subset, methylated==1), .tidysc))

}


# helper fn
.tidysc <- function(sc) {

  sc$frag_grp <- .tidyfg(rownames(sc))
  sc$read_count <- sc$mapped  
  sc[, c("frag_grp", "id", "read_count")]

}


# helper fn
.tidyfg <- function(frags) {

  gsub("[a-z]", "", tolower(sapply(strsplit(frags, "\\-"), `[`, 1)))
  
}
