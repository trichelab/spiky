#' spike-in counts for two samples, as a wide data.frame
#'
#' A data.frame with spike-in results from control samples in the manuscript.
#' This maps 1:1 onto `spike_read_counts` using reshape2::melt.
#'
#' @importFrom utils data
#' @usage data(dedup)
#' @format A data.frame object with
#' \describe{
#'   \item{frag_grp}{the encoded spike contig name: basepairs_CpGs_GCpercent}
#'   \item{read_count_6547}{read coverage for this spike in sample 6547}
#'   \item{read_count_6548}{read coverage for this spike in sample 6548}
#' }
#' @source This data was created using inst/script/loadDedup.R
"dedup"


#' various mitochondrial genomes sometimes used as endogenous spike-ins
#'
#' A DataFrame with species, genome, accession, and sequence for GenBank
#' mitochondrial genome depositions. No concentration provided; add if needed.
#'
#' @usage data(genbank_mito)
#' @format A DataFrame object with
#' \describe{
#'   \item{species}{the species whence the record came, as a character string}
#'   \item{genome}{the genome assembly whence the mtDNA, as a character string}
#'   \item{accession}{the genbank accession, as a character string}
#'   \item{sequence}{genome sequence, as a DNAStringSet}
#' }
#' @source www.ncbi.nlm.nih.gov/genbank/
"genbank_mito"


#' lambda and phiX phage sequences, sometimes used as spike-ins
#'
#' A DataFrame with sequence, methylated, CpGs, GCfrac, and OECpG for phages
#'
#' @usage data(phage)
#' @format A DataFrame object with
#' \describe{
#'   \item{sequence}{genome sequence, as a DNAStringSet}
#'   \item{methylated}{whether CpGs are methylated, as an integer}
#'   \item{CpGs}{the number of CpGs in the phage genome, as an integer}
#'   \item{GCfrac}{the GC fraction of the phage genome, as a numeric}
#'   \item{OECpG}{the observed / expected CpG fraction, as a numeric}
#' }
#' @source www.ncbi.nlm.nih.gov/genbank/
"phage"


#' spike-in contig properties for Sam's cfMeDIP spikes
#'
#' A DataFrame with sequence, concentration, and other properties of Sam's
#' synthetic cfMeDIP spike-in controls. The row names redudantly encode some
#' of these properties, such as the number of CpGs in the spike-in sequence.
#'
#' @usage data(spike)
#' @format A DataFrame object with
#' \describe{
#'   \item{sequence}{contig sequence, as a DNAStringSet}
#'   \item{methylated}{are the CpGs in this spike-in methylated? 0 or 1}
#'   \item{CpGs}{number of CpG dinucleotides in the spike, from 1 to 16}
#'   \item{fmol}{femtomolar concentration of the spike-in for standard mix}
#'   \item{molmass}{molar mass of spike-in sequence}
#' }
#' @source https://doi.org/10.1101/2021.02.12.430289
"spike"


#' spike-in counts, as a long data.frame
#'
#' A data.frame with spike-in results from control samples in the manuscript.
#' This maps 1:1 onto `dedup` using reshape2::melt.
#'
#' @usage data(spike_read_counts)
#' @format A data.frame object with
#' \describe{
#'   \item{frag_grp}{the encoded spike contig name: basepairs_CpGs_GCpercent}
#'   \item{id}{subject from whom cfMeDIP spike reads (column 3) were counted}
#'   \item{read_count}{read coverage for this spike in this subject (column 2)}
#' }
#' @source This data was created using inst/script/loadDedup.R
"spike_read_counts"

#' spike-in counts, as a long data.frame
#'
#' A data.frame with spike-in results from CRAM files (generated from scan_spike_counts(CRAMs, spike=spike))
#'
#' @usage data(spike_cram_counts)
#' @format A data.frame object with
#' \describe{
#'   \item{frag_grp}{the encoded spike contig name: basepairs_CpGs_GCpercent}
#'   \item{id}{subject from whom cfMeDIP spike reads (column 3) were counted}
#'   \item{read_count}{read coverage for this spike in this subject (column 2)}
#' }
#' @source Generated from scan_spike_counts(CRAMs, spike=spike) using example CRAMs
#' containing spike contigs
"spike_cram_counts"


#' scan_spiked_bam results from a merged cfMeDIP CRAM file (chr22 and spikes)
#'
#' A CompressedGRangesList object with `genomic` (chr22) and `spikes` coverage,
#' binned every 300bp for the genomic contigs then averaged across the bin, and
#' summarized for each spike contig as (the default) `max` coverage. (In other
#' words, the default output of scan_spiked_bam, restricted to a small enough
#' set of genomic regions to be practical for examples.) This represents what
#' most users will want to generate from their own merged BAMs or CRAMs, and
#' is used repeatedly in downstream examples throughout the package.
#'
#' @usage data(ssb_res)
#' @format A CompressedGRangesList of coverage results, containing
#' \describe{
#'   \item{genomic}{a GRanges with one metadata column, `coverage`}
#'   \item{spikes}{a GRanges with one metadata column, `coverage`}
#' }
#' @source Generated using scan_spiked_bam on an example bam
#' containing chr22 and spike contigs.
"ssb_res"

#' A Granges object with genomic coverage from chr21q22,
#' binned every 300bp for the genomic contigs then averaged across the bin. (In other
#' words, the default output of scan_genomic_contigs or scan_genomic_bedpe,
#' restricted to a small enough set of genomic regions to be practical for examples.)
#' This represents what most users will want to generate from their own genomic BAMs or BEDPEs, and
#' is used repeatedly in downstream examples throughout the package.
#'
#' @usage data(genomic_res)
#' @format A GRanges of coverage results with one metadata column, `coverage`
#' @source Generated using scan_genomic_bedpe or scan_genomic_contigs on an example bedpe or bam
#' containing chr21q22 contigs.
"genomic_res"

#' A Granges object with spike-in sequence coverage, and
#' summarized for each spike contig as (the default) `max` coverage. (In other
#' words, the default output of scan_spike_contigs or scan_spike_bedpe)
#' This represents what most users will want to generate from their own spike-in BAMs or BEDPEs, and
#' is used repeatedly in downstream examples throughout the package.
#'
#' @usage data(spike_res)
#' @format A GRanges of coverage results with one metadata column, `coverage`
#' @source Generated using scan_spike_bedpe or scan_spike_contigs on an example bedpe or bam
#' containing spike contigs.
"spike_res"

#' a test GRanges with UMI'ed genomic sequences used as controls
#'
#' Sources and overlap widths of various read sequences in a test CRAM.
#'
#' @usage data(testGR)
#' @format A GRanges object with an mcols() DataFrame containing
#' \describe{
#'   \item{UMI1}{the unique molecular identifier on the forward read}
#'   \item{UMI2}{the unique molecular identifier on the reverse read}
#'   \item{seq}{the sequence of the fragment}
#'   \item{name}{the name of the fragment}
#'   \item{score}{whether the fragment passes filters (always 1)}
#' }
#' @source Generated using inst/script/loadTest.R
"testGR"
