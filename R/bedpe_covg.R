#' get coverage from a BEDPE file (possibly removing duplicated bits)
#' 
#' @param x         a Tabixed BEDPE file, or a TabixFile of one
#' @param param     a GRanges of regions to read, or NULL (NULL)
#' @param stranded  treat pairs as stranded? (FALSE)
#' @param undup     remove overlapping bits? (TRUE)
#'
#' @return          an RleList of coverage
#'
#' @import          S4Vectors
#'
#' @export
bedpe_covg <- function(x, param=NULL, stranded=FALSE, undup=TRUE) { 

  p <- read_bedpe(x, param=param, stranded=stranded)
  if (undup) { 
    stop("Need to reconsider this") 
  }
  stop("Finish coverage on BEDPE/BED/fragment files!")
   
}


# helper
.undup <- function(p) { 


}
