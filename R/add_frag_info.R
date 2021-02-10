#' decode fragment identifiers for spike-in standards
#' 
#' given a vector of fragment identifiers like `160_2_35`, encoded as
#' lengthInBp_numberOfCpGs_GCpercent, pull them apart into a data.frame,
#' add those columns to the source data.frame (whence `frag_grp` or similar), 
#' and return the updated data.frame. 
#' 
#' @param x         data.frame with a column of frag_grp information (see above)
#' @param frag_grp  the column name for the frag_grp information (`frag_grp`)
#' 
#' @return          the data.frame x, augmented with fraglen, CpG, and GC cols
#' 
#' @examples 
#' data(spike_read_counts) 
#' add_frag_info(spike_read_counts)
#'
#' @export 
add_frag_info <- function(x, frag_grp="frag_grp") { 

  fg <- gsub("[a-z]+", "", x[, frag_grp], ignore.case=TRUE)
  frag_grp_parse <- data.frame(do.call(rbind, strsplit(fg, "_")))
  names(frag_grp_parse) <- c("fraglen", "CpG", "GC")
  for (i in names(frag_grp_parse)) x[, i] <- as.integer(frag_grp_parse[, i])
  return(x) 

}
