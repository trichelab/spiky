#' take an Rle[List] and either leave it alone (ret="rle") or transform it    
#' 
#' Since a lot of spiky functions have to summarize coverage, this shims that.
#' 
#' @param x         the BAM or CRAM (or BEDPE)
#' @param contigs   the contigs to format
#' @param param     a ScanBamParam object (ignored if BEDPE)
#' @param ret       RleList ("rle") or data.frame ("df")
#' 
#' @return          depends on the value of `ret`
#' 
#' @export 
format_coverage <- function(x, contigs, param=NULL, ret=c("rle", "df")) { 

  bf <- try(BamFile(x))
  ret <- match.arg(ret)

  if (ret == "df") { 
    # return a data.frame
  } else {
    # return the raw Rle
  }
  # format_coverage(bam=bam, param=param, contigs=contigs, ret=ret)
  covg <- coverage(BamFile(bam), param=param)[contigs]


}
