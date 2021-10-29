#' read a BEDPE file into Pairs of GRanges (as if a GAlignmentPairs or similar)
#' 
#' @param x         a Tabixed BEDPE file, or a TabixFile of one
#' @param ...       additional arguments to pass to scanTabix internally
#' @param stranded  is this BEDPE file stranded? (FALSE) 
#' @param fraglen   compute the fragment length? (TRUE) 
#' 
#' @return          a Pairs of GRanges, perhaps with $score or $fraglen
#'
#' @details         BEDPE import in R is a shambles. This is a bandaid on a GSW.
#'                  if all(score(pe) == 1), mcols(pe)$score will be dropped.
#' 
#' @examples
#' \dontrun{
#'   bedpe <- "GSM5067076_2020_A64_bedpe.bed.gz"
#'   WT1_hg38 <- GRanges("chr11", IRanges(32387775, 32435564), "-")
#'   read_bedpe(bedpe, param=WT1_hg38)
#'   # note that since all(score(bedpe) == 1, score is dropped)
#' }
#' 
#' @seealso         bedpe_covg
#'
#' @import          Rsamtools
#' @import          S4Vectors
#'
#' @export
read_bedpe <- function(x, ..., stranded=FALSE, fraglen=TRUE) { 

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

  # drop strand fields if the BEDPE is not stranded 
  if (!stranded) fmt <- fmt[!grepl("strand[12]", names(fmt))]

  # switch to scanTabix -> parsing 
  bits <- list(
    first=fmt[grep(1, names(fmt))],
    second=fmt[grep(2, names(fmt))],
    score=fmt[grep("score", names(fmt))]
  )

  # switch to scanTabix? this is slow AF as implemented
  tbx <- unlist(scanTabix(x, ...), use.names=FALSE)
  res <- lapply(bits, .parse_tabix, fmt=fmt, tbx=tbx)

  # paired 
  pe <- with(res,
             Pairs(first=first, 
                   second=second, 
                   score=as.numeric(score)))
  if (fraglen) mcols(pe)$fraglen <- (end(second(pe)) - start(first(pe)))
  if (all(score(pe) == 1)) mcols(pe)$score <- NULL 
  return(pe) 

}


# helper
.parse_tabix <- function(bits, fmt, tbx) {

  # field positions
  pos <- vapply(names(bits), grep, x=names(fmt), integer(1))
  dat <- do.call(rbind, strsplit(tbx, "\t"))[, pos]

  if (length(bits) > 1) { 
    colnames(dat) <- sub("[12]$", "", names(bits))
    if (all(any(grepl("seqnames", names(bits))),
            any(grepl("start", names(bits))),
            any(grepl("end", names(bits))))) {
      dat <- makeGRangesFromDataFrame(data.frame(dat))
    }
  }

  return(dat)

}


# helper
.rename <- function(x, suffix="[12]$") {

  names(x) <- sub(suffix, "", names(x))
  return(x) 

}
