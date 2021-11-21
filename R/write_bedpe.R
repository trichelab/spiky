#' write out a Pairs object as a BEDPE file
#' 
#' see read_bedpe for more details on the BEDPE format itself; note that this 
#' could eventually be patched into rtracklayer, but don't hold your breath. 
#'
#' @param object  a Pairs object (contrast with most export... functions)
#' @param con     a connection (currently just a filename, see above)
#' 
#' @return        the results of write.table(..., con)
#' 
#' @details       Unlike read_bedpe, write_bedpe tries to respect extra columns
#' 
#' @import        GenomicRanges
#' 
#' @export
write_bedpe <- function(object, con, ...) { 
 
  # use S4 generics to fix 
  if (!is(object, "Pairs")) {
    stop("write_bedpe currently requires a Pairs object as its first argument.")
  }
  
  # use S4 generics to fix 
  if (!is(con, "character") | (length(con) > 1)) {
    stop("write_bedpe currently takes a single string as its second argument.")
  }

  # disgusting yet effective 
  grcols <- c("seqnames", "start", "end")
  firstcols <- paste0("first.", grcols)
  secondcols <- paste0("second.", grcols) 
  required <- c(firstcols, secondcols)

  # eventually this should happen to read_bedpe, too
  extras <- c("name", "score", "first.strand", "second.strand")
  optional <- intersect(extras, names(mcols(object)))
  optional <- c(optional, setdiff(names(mcols(object)), c(required, extras)))

  # cf. Mitch Guttman's SPRITE files
  reordered <- c(required, optional) 

  # high-tech output function
  # maybe it would be worth sorting and tabixing this
  write.table(as.data.frame(object)[, reordered], file=con,
              row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

}
