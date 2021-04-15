#' tabulate methylation specificity for multiple spike-in BAM/CRAM files
#'
#' Methylation specificity is here defined as methylated_spike_covg/spike_covg
#'
#' @param files   a vector of BAM/CRAM file names
#' @param spike   a spike-in database
#' @param sep     the separator for spike-in contig names ("_")
#'
#' @return        a matrix with columns "mean" and "median"
#'
#' @examples
#' data(spike)
#' library(GenomicRanges)
#' sb <- system.file("extdata", "example.spike.bam", package="spiky",
#'                   mustWork=TRUE)
#' scan_methylation_specificity(sb, spike=spike)
#'
#' @export
scan_methylation_specificity <- function(files, spike=NULL, sep="_") {
  if (is.null(spike)) spike = spiky::spike
  stopifnot(all(file.exists(files)))
  ssc <- lapply(files, spike_counts, spike=spike, dump_idx=TRUE)
  if (is.null(names(ssc))) {
    names(ssc) <- gsub("\\.(bam|cram)", "", basename(files))
  }
  t(vapply(ssc, .tally, numeric(2)))

}


# helper fn
.tally <- function(sc) {

  methylated <- rowsum(subset(sc, methylated==1)$mapped,
                       .rfg(subset(sc, methylated==1)))
  total <- rowsum(sc$mapped, .rfg(sc))
  c(mean=mean(methylated/total), median=median(methylated/total))

}


# helper fn (row frag grp)
.rfg <- function(x) {

  gsub("[a-z]", "", tolower(sapply(strsplit(rownames(x), "\\-"), `[`, 1)))

}
