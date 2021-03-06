% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_spikes.R
\name{process_spikes}
\alias{process_spikes}
\title{QC, QA, and processing for a new spike database}
\usage{
process_spikes(fasta, methylated = 0, ...)
}
\arguments{
\item{fasta}{fasta file (or GRanges or DataFrame) w/spike sequences}

\item{methylated}{whether CpGs in each are methylated (0 or 1, default 0)}

\item{...}{additional arguments, e.g. kernels (currently unused)}
}
\value{
\preformatted{         a DataFrame suitable for downstream processing
}
}
\description{
Sequence feature verification: never trust anyone, least of all yourself.
}
\details{
GCfrac is the GC content of spikes as a proportion instead of a percent.
OECpG is (observed/expected) CpGs (expectation is 25\% of GC dinucleotides).
}
\examples{

data(spike)
spikes <- system.file("extdata", "spikes.fa", package="spiky", mustWork=TRUE)
spikemeth <- spike$methylated
process_spikes(spikes, spikemeth)

data(phage)
phages <- system.file("extdata", "phages.fa", package="spiky", mustWork=TRUE)
identical(process_spikes(phage), phage)
identical(phage, process_spikes(phage))

data(genbank_mito)
(mt <- process_spikes(genbank_mito)) # see also genbank_mito.R
gb_mito <- system.file("extdata", "genbank_mito.R", package="spiky")

\donttest{
  library(kebabs)
  CpGmotifs <- c("[AT]CG[AT]","C[ATC]G", "CCGG", "CGCG")
  mot <- motifKernel(CpGmotifs, normalized=FALSE)
  km <- getKernelMatrix(mot, subset(phage, methylated == 0)$sequence)
  heatmap(km, symm=TRUE)

  #'   # refactor this out
  mt <- process_spikes(genbank_mito)
  mtiles <- unlist(tileGenome(seqlengths(mt$sequence), tilewidth=100))
  bymito <- split(mtiles, seqnames(mtiles))
  binseqs <- getSeq(mt$sequence, bymito[["Homo sapiens"]])
  rCRS6mers <- kmers(binseqs, k=6)

  # plot binned Pr(kmer):
  library(ComplexHeatmap)
  Heatmap(kmax(rCRS6mers), name="Pr(kmer)")

  # not run
  library(kebabs)
  kernels <- list(
    k6f=spectrumKernel(k=6, r=1),
    k6r=spectrumKernel(k=6, r=1, revComp=TRUE)
  )

  kms <- lapply(kernels, getKernelMatrix, x=mt["Human", "sequence"])
  library(ComplexHeatmap)
}


}
\seealso{
kmers
}
