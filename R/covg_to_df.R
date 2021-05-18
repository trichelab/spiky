#' reshape `scan_spiked_bam` results into data.frames for model_glm_pmol
#' 
#' @param res   the result of running scan_spiked_bam
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
#' data(ssb_res, package="spiky") 
#' subsetted <- covg_to_df(ssb_res, spike=spike, meth=TRUE) 
#' summed <- covg_to_df(ssb_res, spike=spike, meth=FALSE)
#' round((summed$read_count - subsetted$read_count) / summed$read_count, 3)
#' 
#' @export 
covg_to_df <- function(res, spike, meth=TRUE, ID=NULL) {

  if (is.null(ID)) { # mostly to keep from going insane in testing
    ID <- paste(toupper(sample(c(letters, 0:9), 4)), collapse="")
  }

  if ("spikes" %in% names(res)) {
    spikes <- res$spikes
  } else {
    spikes <- res
  }
  if (!"methylated" %in% names(mcols(spikes))) {
    if (is.null(spike)) {
      message("You need to provide a `spike` database to proceed.")
      message("Consider using the defaults, i.e. `data(spike, package='spiky')`.")
      stop("Cannot assign methylated/unmethylated status to spike contigs.")
    }
    spikes$methylated <- spike[names(spikes), "methylated"]
  }
  spikes$stub <- vapply(strsplit(names(spikes), "\\-"), `[`, character(1), 1)
  spikes <- as.data.frame(mcols(spikes[, c("stub","coverage","methylated")]))
  spikes$id <- ID 

  if (meth) {
    spikes <- subset(spikes, spikes$methylated == 1)
    spikes <- spikes[, c("stub", "id", "coverage")] 
    names(spikes) <- c("frag_grp", "id", "read_count")
  } else {
    spikes <- data.frame(id=ID, 
                         read_count=tapply(spikes$coverage, spikes$stub, sum))
    spikes$frag_grp <- rownames(spikes) 
    spikes <- spikes[, c("frag_grp", "id", "read_count")] 
  }

  return(spikes[order(spikes$frag_grp), c("frag_grp", "id", "read_count")])

}
