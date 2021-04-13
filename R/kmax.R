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
#' # not run
#' # library(BSgenome.Mmusculus.UCSC.mm10.masked)
#' # mm10k6 <- kmers(Mmusculus)
#' # rownames(mm10k6) <- paste0("mm10_", rownames(mm10k6))
#' # 
#' # library(BSgenome.Hsapiens.UCSC.hg38.masked)
#' # hg38k6 <- kmers(Hsapiens) 
#' # rownames(hg38k6) <- paste0("hg38_", rownames(hg38k6))
#' #
#' # hgmmphmtk6 <- rbind(hg38k6[paste0("hg38_chr", 1:22), ],
#' #                     mm10k6[paste0("mm10_chr", 1:19), ],
#' #                     phk6, mtk6)
#' #
#' # library(ComplexHeatmap)
#' # Heatmap(kmax(hgmmphmtk6), name="Pr(kmer)")
#' 
#' @export 
kmax <- function(km, normalize=TRUE) {

  if (length(unique(nchar(colnames(km)))) > 1) stop("kmax is for kmer tallies!")
  if (normalize) km <- sweep(km, 1, rowSums(km), `/`)
  km[seq_len(nrow(km)), unique(max.col(km))]

}
