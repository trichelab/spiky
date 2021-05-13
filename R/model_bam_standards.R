#' Build a Bayesian additive model from spike-ins to correct bias in *-seq
#'
#' @param   x     data with assorted feature information (GCfrac, CpGs, etc)
#' @param   conc  concentration for each spike (must be provided!)
#' @param   fm    model formula (conc ~ read_count + fraglen + GCfrac + CpGs_3)
#' @param   ...   other arguments to pass to `bamlss`
#'
#' @return        the model fit for the data
#'
#' @examples
#'
#' library(bamlss)
#' data(spike_cram_counts,package="spiky")
#' data(spike,package="spiky")
#' scc <- add_frag_info(spike_cram_counts, spike=spike)
#' scc$conc <- scc$conc * 0.9 # adjust for dilution
#' scc$CpGs_3 <- scc$CpGs ^ (1/3)
#' fit0 <- model_bam_standards(scc,
#'                             fm=conc ~ read_count + fraglen)
#' fit1 <- model_bam_standards(scc,
#'                             fm=conc ~ read_count + fraglen + GCfrac + CpGs_3)
#' DIC(fit0, fit1)
#'
#' @import        bamlss
#'
#' @export
model_bam_standards <- function(x, conc=NULL, fm=NULL, ...) {

  if (is.null(fm)) fm <- conc ~ read_count + fraglen + GCfrac + CpGs_3
  if (is.null(x$conc) & !is.null(conc)) x$conc <- conc

  cols <- strsplit(as.character(fm)[3], " + ", fixed=TRUE)[[1]]
  if (!all(cols %in% names(x))) {
    stop("Input data must have columns ", paste(cols, collapse=", "))
  }

  # Gaussian model by default -- can alter with ... params
  bamlss(fm, data=x, ...)

}
