#' QC, QA, and processing for a new spike database
#'
#' Sequence feature verification: never trust anyone, least of all yourself.
#'
#' @param fasta         fasta file (or GRanges or DataFrame) w/spike sequences
#' @param methylated    whether CpGs in each are methylated (0 or 1, default 0)
#' @param ...           additional arguments, e.g. kernels (currently unused)
#'
#' @return              a DataFrame suitable for downstream processing
#'
#' @details
#' GCfrac is the GC content of spikes as a proportion instead of a percent.
#' OECpG is (observed/expected) CpGs (expectation is 25% of GC dinucleotides).
#'
#' @examples
#'
#' data(spike)
#' spikes <- system.file("extdata", "spikes.fa", package="spiky", mustWork=TRUE)
#' spikemeth <- spike$methylated
#' process_spikes(spikes, spikemeth)
#'
#' data(phage)
#' phages <- system.file("extdata", "phages.fa", package="spiky", mustWork=TRUE)
#' identical(process_spikes(phage), phage)
#' identical(phage, process_spikes(phage))
#'
#' data(genbank_mito)
#' (mt <- process_spikes(genbank_mito)) # see also genbank_mito.R
#' gb_mito <- system.file("extdata", "genbank_mito.R", package="spiky")
#'
#' \donttest{
#'   library(kebabs)
#'   CpGmotifs <- c("[AT]CG[AT]","C[ATC]G", "CCGG", "CGCG")
#'   mot <- motifKernel(CpGmotifs, normalized=FALSE)
#'   km <- getKernelMatrix(mot, subset(phage, methylated == 0)$sequence)
#'   heatmap(km, symm=TRUE)
#'
#'   #'   # refactor this out
#'   mt <- process_spikes(genbank_mito)
#'   mtiles <- unlist(tileGenome(seqlengths(mt$sequence), tilewidth=100))
#'   bymito <- split(mtiles, seqnames(mtiles))
#'   binseqs <- getSeq(mt$sequence, bymito[["Homo sapiens"]])
#'   rCRS6mers <- kmers(binseqs, k=6)
#'
#'   # plot binned Pr(kmer):
#'   library(ComplexHeatmap)
#'   Heatmap(kmax(rCRS6mers), name="Pr(kmer)")
#'
#'   # not run
#'   library(kebabs)
#'   kernels <- list(
#'     k6f=spectrumKernel(k=6, r=1),
#'     k6r=spectrumKernel(k=6, r=1, revComp=TRUE)
#'   )
#'
#'   kms <- lapply(kernels, getKernelMatrix, x=mt["Human", "sequence"])
#'   library(ComplexHeatmap)
#' }
#'
#'
#' @import Biostrings
#' @import S4Vectors
#'
#' @seealso kmers
#'
#' @export
process_spikes <- function(fasta, methylated=0, ...) {

  if (is.character(fasta) && file.exists(fasta)) {
    spikes <- DataFrame(sequence=readDNAStringSet(fasta), methylated=methylated)
  } else if (is(fasta, "GRanges") && "sequence" %in% names(mcols(fasta))) {
    spikes <- fasta
    if (!"methylated" %in% names(mcols(fasta))) spikes$methylated <- methylated
  } else if (is(fasta, "DataFrame") && "sequence" %in% names(fasta)) {
    spikes <- fasta
    if (!"methylated" %in% names(fasta)) spikes$methylated <- methylated
  } else {
    stop(fasta, " does not exist or is not in a suitable format to proceed.")
  }

  spikes$CpGs <- dinucleotideFrequency(spikes$sequence)[, "CG"]
  GC <- apply(alphabetFrequency(spikes$sequence)[, c("C", "G")], 1, sum)
  spikes$GCfrac <- round(GC / width(spikes$sequence), 2)
  spikes$OECpG <- round(spikes$CpGs / (GC/8), 2) # CpG expected 25% of GC dinuc
  return(spikes)

}
