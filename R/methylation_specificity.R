#' compute methylation specificity for spike-in standards
#' 
#' In a cfMeDIP experiment, the yield of methylated fragments should be >95% 
#' (ideally 98-99%) due to the nature of the assay. 
#'
#' @export
methylation_specificity <- function(spike_data, ...) { 
  data("spike")
  methylated_spikes <- rownames(spike[spike$methylated == 1,])
  sum(spike_data[methylated_spikes])
  methyl_spec <- sum(spike_data[methylated_spikes]) / (sum(spike_data))
  return(methyl_spec)
}
