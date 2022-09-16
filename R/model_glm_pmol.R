#' Build a generalized linear model from spike-ins to correct bias in cfMeDIP
#'
#' formerly '2020_model_glm_fmol'.  Note that everything in x can be had from
#' a BAM/CRAM with spike contigs named as frag_grp (len_CpGs_GC) in the index
#' and in fact that is what scan_spiked_bam now does.
#'
#' @param   x     data w/frag_grp, id, and read_count; or scan_spiked_bam result
#' @param   conc  concentration for each spike (will be referenced if NULL)
#' @param   spike spike database, e.g. `data(spike, package='spiky')`
#' @param   ...   other arguments to pass to `glm` (e.g. `family`)
#'
#' @return        the model fit for the data
#'
#' @examples
#'
#' data(spike, package="spiky")
#'
#' data(spike_read_counts, package="spiky")
#' fit1 <- model_glm_pmol(spike_read_counts, spike=spike)
#'
#' data(spike_res) # scan_spiked_bam result
#' fit2 <- model_glm_pmol(spike_res, spike=spike)
#'
#' @importFrom stats glm
#'
#' @export
model_glm_pmol <- function(x, spike, conc=NULL, ...) {

  if (!is(x, "data.frame")) {
    message("Converting x (a coverage result?) to a data.frame...")
    x <- covg_to_df(x, spike=spike)
  }

  cols <- c("frag_grp", "id", "read_count")
  if (!all(cols %in% names(x))) {
    stop("Input data must have columns frag_grp, id, and read_count")
  }

  x <- add_frag_info(x, "frag_grp") # now exported
  x$conc <- .getConcFromFraglen(x$fraglen,conc=conc, spike=spike) # below

  # Adjust for the 0.01ng dilution
  x$conc <- x$conc * 0.9

  # Cube root CpG number to normalize
  x$CpG <- as.integer(x$CpG) # this too would be helpful to standardize
  x$CpG_3 <- x$CpG ^ (1/3)

  # Gaussian model by default -- can alter with ... params
  fit <- glm(formula = conc ~ read_count + fraglen + GC + CpG_3, data = x, ...)
  r2_gaussian= 1 - (fit$deviance / fit$null.deviance)
  attr(fit, "r2_gaussian") <- r2_gaussian #  report in summary
  attr(fit, "data") <- x
  return(fit)

}


# helper fn
.getConcFromFraglen <- function(fraglen, spike, conc, concs=NULL) {

  # pull these from spike
  if (is.null(concs)) {
    tbl <- with(spike, table(width(sequence), conc))
    concs <- as.numeric(colnames(tbl)[max.col(tbl)])
    names(concs) <- rownames(tbl)
  }

  # show(concs)
  #    80   160   320
  # 0.004 0.002 0.001

  res <- concs[as.character(fraglen)]

  # essentially the fallthrough from ifelse
  res[is.na(res)] <- concs[length(concs)]

  return(res)

}
