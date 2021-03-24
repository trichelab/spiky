#' decode fragment identifiers for spike-in standards
#' 
#' given a vector of fragment identifiers like `160_2_35` or `80b_1C_35G-2`,
#' encoded typically as lengthInBp_numberOfCpGs_GCpercent, and optionally a 
#' database of spike-in sequences corresponding to those fragments, add those
#' columns to the source data (along with, if present in the database, other 
#' metadata such as standard concentrations, GC fraction, etc.) and return i
#' an updated DataFrame. 
#' 
#' @param x         data.frame with a column of spike information (see above)
#' @param frag_grp  column name for the spike contig information (`frag_grp`)
#' @param spike     optional database of spike-in properties (none)
#' 
#' @return          the data.frame x, augmented with metadata columns
#' 
#' @examples 
#' data(spike) 
#' data(spike_cram_counts) 
#' spike <- subset(spike, methylated == 1)
#' add_frag_info(spike_cram_counts, spike=spike)
#'
#' @export 
add_frag_info <- function(x, frag_grp="frag_grp", spike=NULL) { 

  fg <- gsub("[a-z]+", "", x[, frag_grp], ignore.case=TRUE)
  if (!is.null(spike) && (!frag_grp %in% names(spike))) {
    spike[, frag_grp] <- .tidyfg(rownames(spike))
  }
  frag_grp_parse <- data.frame(do.call(rbind, strsplit(fg, "_")))
  names(frag_grp_parse) <- c("fraglen", "CpG", "GC")
  for (i in names(frag_grp_parse)) x[, i] <- as.integer(frag_grp_parse[, i])
  if (!is.null(spike)) { # can override the above
    for (name in setdiff(names(spike), "sequence")) {
      x[, name] <- spike[match(x[, frag_grp], spike[, frag_grp]), name]
    }
  }
  return(x) 

}
