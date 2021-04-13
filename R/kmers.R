#' oligonucleotideFrequency, but less letters and more convenient.
#' 
#' @param  x          BSgenome, DFrame with `sequence` column, or DNAStringSet
#' @param  k          the length of the kmers (default is 6)
#' 
#' @return            a matrix of contigs (rows) by kmer frequencies (columns)
#'
#' @details 
#' The companion `kmax` function finds the maximum frequency kmer for each
#' contig and plots all of them together for comparison purposes. 
#' 
#' @examples
#' 
#' data(genbank_mito, package="spiky") 
#' mtk6 <- kmers(genbank_mito, k=6)
#' kmax(mtk6)
#' 
#' data(phage, package="spiky") 
#' phk6 <- kmers(phage, k=6)
#' kmax(phk6)
#' 
#' @seealso           kmax
#' 
#' @import            Biostrings
#' 
#' @export
kmers <- function(x, k=6) {

  if (is(x, "DNAStringSet")) {
    res <- oligonucleotideFrequency(x, width=k)
    rownames(res) <- names(x)
  } else if (is(x, "DFrame")) { 
    res <- oligonucleotideFrequency(x$sequence, width=k)
    rownames(res) <- rownames(x)
  } else if (is(x, "BSgenome")) {
    chrs <- names(x)
    names(chrs) <- chrs
    res <- t(sapply(chrs, 
                    function(chr) oligonucleotideFrequency(x[[chr]], k)))
  } else { 
    stop("Don't know how to get kmers for an object of type ", class(x))
  } 

  attr(res, "what") <- "kmers"
  return(res)

}
