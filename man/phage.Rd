% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{phage}
\alias{phage}
\title{lambda and phiX phage sequences, sometimes used as spike-ins}
\format{
A DataFrame object with
\describe{
\item{sequence}{genome sequence, as a DNAStringSet}
\item{methylated}{whether CpGs are methylated, as an integer}
\item{CpGs}{the number of CpGs in the phage genome, as an integer}
\item{GCfrac}{the GC fraction of the phage genome, as a numeric}
\item{OECpG}{the observed / expected CpG fraction, as a numeric}
}
}
\source{
www.ncbi.nlm.nih.gov/genbank/
}
\usage{
data(phage)
}
\description{
A DataFrame with sequence, methylated, CpGs, GCfrac, and OECpG for phages
}
\keyword{datasets}
