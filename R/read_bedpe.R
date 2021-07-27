#' read a BEDPE file into Pairs of GRanges (as if a GAlignmentPairs or similar)
#' 
#' @param x         a Tabixed BEDPE file, or a TabixFile of one
#' @param param     a GRanges of regions to read, or NULL (NULL)
#' @param stranded  retain strand of GRanges in the resulting Pairs? (TRUE) 
#' 
#' @return          an S4Vectors::Pairs of GRanges with $score in mcols() 
#'
#' @details         BEDPE import in R is a shambles. This is a bandaid on a GSW.
#' 
#' @seealso         bedpe_covg
#'
#' @import          Rsamtools
#' @import          S4Vectors
#'
#' @export
read_bedpe <- function(x, param=NULL, stranded=TRUE) { 

  # mandatory for proper scanning and param support
  if (!is(x, "TabixFile")) x <- TabixFile(x)

  # format specification for BEDPE 
  fmt <- c(seqnames1="character", 
           start1="integer", 
           end1="integer",
           seqnames2="character", 
           start2="integer", 
           end2="integer",
           strand1="character",
           strand2="character",
           score="numeric")

  # switch to scanTabix -> parsing 
  bits <- list(
    first=fmt[grep(1, names(fmt))],
    second=fmt[grep(2, names(fmt))],
    score=fmt[grep("score", names(fmt))]
  )
  bits <- lapply(bits, .rename) # sanitize

  # switch to scanTabix? this is slow AF as implemented
  res <- lapply(bits, 
                .parse_tabix, 
                param=param, 
                stranded=stranded,
                x=scanTabix(x, param=param))

  with(res,
       Pairs(first=first, 
             second=second, 
             score=score))

}


# helper
.parse_tabix <- function(x, fields, param=NULL, stranded=TRUE) {

  if (all(any(grepl("seqnames", names(fields))),
          any(grepl("start", names(fields))),
          any(grepl("end", names(fields))))) {
  
    # a GRanges is being requested
    # drop the strand if !stranded
    if (!stranded) fields <- fields[!grepl("strand", names(fields))]

    n <- length(fields)
    res <- data.frame(t(vapply(strsplit(x, "\t"), `[`, character(n), fields)))
    names(res) <- names(fields)

  } else { 

    # hard to say

  }

  gr <- makeGRangesFromDataFrame(res, keep=TRUE) 
  if (!stranded) gr <- unstrand(gr)
  return(gr)

}


# helper
.rename <- function(x, suffix="[12]$") {

  names(x) <- sub(suffix, "", names(x))
  return(x) 

}
