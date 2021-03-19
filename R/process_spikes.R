#' QC, QA, and processing for a new spike database 
#' 
#' Never trust anyone, least of all yourself 
#'
#' @param fasta         fasta file with new spike sequences 
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
#' phages <- system.file("extdata", "phages.fa", package="spiky", mustWork=TRUE)
#' process_spikes(phages)
#' 
#' data(spike)
#' spikes <- system.file("extdata", "spikes.fa", package="spiky", mustWork=TRUE)
#' spikemeth <- spike$methylated
#' process_spikes(spikes, spikemeth)
#' 
#' # not run
#' # library(kebabs)
#' # CpGmotifs <- c("[AT]CG[AT]","C[ATC]G", "CCGG", "CGCG")
#' # mot <- motifKernel(CpGmotifs, normalized=FALSE)
#' # km <- getKernelMatrix(mot, subset(phages, methylated == 0)$sequence)
#' # heatmap(km, symm=TRUE)
#' 
#' @import Biostrings
#' @import S4Vectors
#' 
#' @export
process_spikes <- function(fasta, methylated=0, ...) {

  stopifnot(file.exists(fasta))
  spikes <- DataFrame(sequence=readDNAStringSet(fasta), methylated=methylated)
  spikes$CpGs <- dinucleotideFrequency(spikes$sequence)[, "CG"]
  GC <- apply(alphabetFrequency(spikes$sequence)[, c("C", "G")], 1, sum)
  spikes$GCfrac <- round(GC / width(spikes$sequence), 2)
  spikes$OECpG <- round(spikes$CpGs / (GC/8), 2) # CpG expected 25% of GC dinuc
  return(spikes)

}
