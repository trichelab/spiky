#' refactored out of rename_spikes and rename_seqlevels
#' 
#' A common task between generate_spike_fasta, rename_spikes, and 
#' rename_seqlevels is to determine what the largest common subset of 
#' characters between existing contig names and stored standardized contigs 
#' might be. This function eases that task.
#' 
#' @param     contig_names    the names of contigs
#' @param     sep             separator character in contig names ("_")
#' 
#' @return                    a vector of atoms 1:3 from each contig name 
#' 
#' @export
get_base_name <- function(contig_names, sep="_") {

  sapply(lapply(strsplit(contig_names, sep), `[`, 1:3), paste, collapse=sep)

}
