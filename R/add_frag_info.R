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
#' @export 
add_frag_info <- function(x, frag_grp="frag_grp") { 

  frag_grp_parse <- data.frame(do.call(rbind, strsplit(x[, frag_grp], "_")))
  names(frag_grp_parse) <- c("fraglen", "CpG", "GC")
  for (i in names(frag_grp_parse)) x[, i] <- as.integer(frag_grp_parse[, i])
  return(x) 

}


# helper function
.getConcFromFraglen <- function(fraglen, concs = NULL) {
  
  if (is.null(concs)) concs <- c("80" = 0.004, "160" = 0.002, "320" = 0.001)
  res <- concs[as.character(fraglen)]

  # essentially the fallthrough from ifelse
  res[is.na(res)] <- concs[length(concs)]

  return(res)

}
