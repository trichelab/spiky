#' simple contig kmer comparisons
#'
#' @param km          kmer summary
#' @param normalize   normalize (divide by row sums)? (TRUE)
#'
#' @return            the most common kmers for each contig, across all contigs
#'
#' @examples
#'
#'
#' data(genbank_mito, package="spiky")
#' mtk6 <- kmers(genbank_mito, k=6)
#' rownames(mtk6) <- paste0(rownames(mtk6), "_MT")
#' kmax(mtk6)
#'
#' data(phage, package="spiky")
#' phk6 <- kmers(phage, k=6)
#' kmax(phk6, normalize=FALSE)
#'
#' stopifnot(identical(colnames(phk6), colnames(mtk6)))
#' k6 <- rbind(mtk6, phk6)
#' kmax(k6)
#'
#' @export
kmax <- function(km, normalize=TRUE) {

  if (length(unique(nchar(colnames(km)))) > 1) stop("kmax is for kmer tallies!")
  if (normalize) km <- sweep(km, 1, rowSums(km), `/`)
  km[seq_len(nrow(km)), unique(max.col(km))]

}
