#' read a BEDPE file into Pairs of GRanges (as if a GAlignmentPairs or similar)
#'
#' @param x         a Tabixed BEDPE file, or a TabixFile of one
#' @param ...       additional arguments to pass to scanTabix internally
#' @param optional  scan the optional columns (name, score, strand1)? (FALSE)
#' @param fraglen   compute the fragment length? (TRUE)
#' @param keep      keep additional columns? (FALSE)
#' @param stranded  Is the data stranded? (FALSE)
#'
#' @return          a Pairs of GRanges, perhaps with $score or $fraglen
#'
#' @details         BEDPE import in R is a shambles. This is a bandaid on a GSW.
#'                  See the \href{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format}{BEDPE format definition} for full details.
#'                  In short, for a pair of ranges 1 and 2, we have fields
#'                  chrom1, start1, end1, chrom2, start2, end2, and (optionally)
#'                  name, score, strand1, strand2, plus any other user defined
#'                  fields that may be included (these are not yet supported
#'                  by read_bedpe). For example, two valid BEDPE lines are:
#'
#'                  chr1  100   200   chr5  5000  5100  bedpe_example1  30
#'                  chr9  900  5000   chr9  3000  3800  bedpe_example2  99  +  -
#'
#' @examples
#' \dontrun{
#'   bedpe <- "GSM5067076_2020_A64_bedpe.bed.gz"
#'   WT1_hg38 <- GRanges("chr11", IRanges(32387775, 32435564), "-")
#'   read_bedpe(bedpe, param=WT1_hg38)
#' }
#'
#' @seealso         bedpe_covg
#'
#' @import          Rsamtools
#' @import          S4Vectors
#'
#' @export
read_bedpe <- function(x, ..., stranded=FALSE, fraglen=TRUE, optional=FALSE, keep=FALSE) {

  # no support yet for additional mcols, even though it's kind of easy
  if (keep) message("Extra columns are not yet supported.")

  # mandatory for proper scanning and param support
  if (!is(x, "TabixFile")) x <- TabixFile(x)

  # format specification for BEDPE
  fmt <- c(seqnames1="character",
           start1="integer",
           end1="integer",
           seqnames2="character",
           start2="integer",
           end2="integer",
           name="name",
           score="numeric",
           strand1="character",
           strand2="character")

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
  if (length(strsplit(tbx[1], "\t")[[1]]) < length(fmt)) { fmt <- fmt[!names(fmt)=="name"]}
  if (length(strsplit(tbx[1], "\t")[[1]]) < length(fmt)) {stop("Missing columns in input file.")}
  res <- lapply(bits, .parse_tabix, fmt=fmt, tbx=tbx)
  pe <- .toPairs(res)

  if (fraglen) mcols(pe)$fraglen <- (end(second(pe)) - start(first(pe)))
  if (all(is.na(mcols(pe)$score)) | all(mcols(pe)$score == 1)) mcols(pe)$score <- NULL

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


# bit of a hack utility
.toGRs <- function(x) {

  cols <- list(first=grep("(chrom|seqnames|start|end)1", colnames(x), val=TRUE),
               second=grep("(chrom|seqnames|start|end)2", colnames(x),val=TRUE))
  lapply(cols, function(y) makeGRangesFromDataFrame(.rename(x[, y])))

}


# helper
.toPairs <- function(res) {

  paired <- with(res, Pairs(first=first, second=second))
  for (extracol in setdiff(names(res), c("first", "second"))) {
    if (extracol == "score") mcols(paired)$score <- as.numeric(res$score)
    else mcols(paired)[, extracol] <- res[[extracol]]
  }
  return(paired)

}
