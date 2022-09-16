#' reshape `scan_spiked_bam` results into data.frames for model_glm_pmol
#' @param   spike_gr    GRanges of spike contigs (e.g. output object from scan_spiked_bam, scan_spike_contigs, or scan_spike_bedpe)
#' @param spike spike database (as from data(spike, package="spiky"))
#' @param meth  only keep methylated spike reads? (TRUE; if FALSE, sum both)
#' @param ID    an identifier for this sample, if running several (autogenerate)
#'
#' @return      a data.frame with columns 'frag_grp', 'id', and 'read_count'
#'
#' @seealso     scan_spiked_bam
#'
#' @examples
#' data(spike, package="spiky")
#' data(spike_res, package="spiky")
#' subsetted <- covg_to_df(spike_res, spike=spike, meth=TRUE)
#' summed <- covg_to_df(spike_res, spike=spike, meth=FALSE)
#' round((summed$read_count - subsetted$read_count) / summed$read_count, 3)
#'
#' @export
covg_to_df <- function(spike_gr, spike, meth=TRUE, ID=NULL) {

  if (is.null(ID)) { # mostly to keep from going insane in testing
    ID <- paste(toupper(sample(c(letters, 0:9), 4)), collapse="")
  }

  if (!"methylated" %in% names(mcols(spike_gr))) {
    if (is.null(spike)) {
      message("You need to provide a `spike` database to proceed.")
      message("Consider using the defaults, i.e. `data(spike, package='spiky')`.")
      stop("Cannot assign methylated/unmethylated status to spike contigs.")
    }
    spike_gr$methylated <- spike[names(spike_gr), "methylated"]
  }
  spike_gr$stub <- vapply(strsplit(names(spike_gr), "\\-"), `[`, character(1), 1)
  spike_gr <- as.data.frame(mcols(spike_gr[, c("stub","coverage","methylated")]))
  spike_gr$id <- ID

  if (meth) {
    spike_gr <- subset(spike_gr, spike_gr$methylated == 1)
    spike_gr <- spike_gr[, c("stub", "id", "coverage")]
    names(spike_gr) <- c("frag_grp", "id", "read_count")
  } else {
    spike_gr <- data.frame(id=ID,
                         read_count=tapply(spike_gr$coverage, spike_gr$stub, sum))
    spike_gr$frag_grp <- rownames(spike_gr)
    spike_gr <- spike_gr[, c("frag_grp", "id", "read_count")]
  }

  return(spike_gr[order(spike_gr$frag_grp), c("frag_grp", "id", "read_count")])

}
