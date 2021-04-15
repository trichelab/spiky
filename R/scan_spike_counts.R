#' run spike_counts on BAM/CRAM files and shape the results for model_glm_pmol
#'
#' Typically one will want to fit a correction model to multiple samples.
#' This function eases this task by merging the output of spike_counts into
#' a data.frame that model_glm_pmol can directly fit.
#'
#' @param files   a vector of BAM/CRAM file names
#' @param spike   a spike-in database
#' @param sep     the separator for spike-in contig names ("_")
#' @param methylated  a logical (1/0) if only methylated sequences should be returned
#'
#' @return        a data.frame with columns "frag_grp", "id", and "read_count"
#'
#' @examples
#' data(spike)
#' library(GenomicRanges)
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' scan_spike_counts(sb, spike=spike,methylated=1)
#' fit <- model_glm_pmol(scan_spike_counts(sb, spike=spike),spike=spike)
#'
#' @export
scan_spike_counts <- function(files, spike=NULL, sep="_",methylated=1) {

  if (is.null(spike)) spike = spiky::spike
  stopifnot(all(file.exists(files)))
  ssc <- lapply(files, spike_counts, spike=spike, dump_idx=TRUE)
  do.call(rbind, lapply(lapply(ssc, subset, methylated=methylated), .tidysc))

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
