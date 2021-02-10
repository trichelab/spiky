#' predict picomoles of DNA from a fit and read counts (coverage)
#' 
#' FIXME: this could be made MUCH faster by precomputing CpG/GC stats per bin
#' 
#' @param   fit       result of model_glm_pmol 
#' @param   ssb_res   the data / new data 
#' @param   bsgenome  BSgenome name (if null, will guess from ssb_res)
#' @param   ret       return a data.frame ("df") or GRanges ("gr")?  ("gr")
#' 
#' @details
#' Using GRanges as the return value is (perhaps counterintuitively) *much* 
#' faster than the data.frame, since the sequence of the bins gets converted
#' from a BSgenome representation to characters in the latter (it is implied
#' by the bin start, stop, and genome when left as a GRanges). 
#' 
#' @return  object with read count, fraglen, GC%, CpG**(1/3), and concentration
#' 
#' @import tools
#' @import BSgenome 
#' @import S4Vectors
#' @import Biostrings
#' @importFrom stats predict
#' 
#' @examples 
#' 
#' data(ssb_res)
#' fit <- model_glm_pmol(covg_to_df(ssb_res)) 
#' preddf <- predict_pmol(fit, ssb_res, ret="df") 
#' predgr <- predict_pmol(fit, ssb_res, ret="gr") 
#' 
#' @export 
predict_pmol <- function(fit, ssb_res, bsgenome=NULL, ret=c("df", "gr")) {

  # what shall be returned?  (eventually it would be nice to pass GRs around)
  ret <- match.arg(ret)

  # guess the genome if not provided (it won't be)
  if (is.null(bsgenome)) {
    assembly <- grep("(hg|GRCh)", unique(genome(ssb_res)), value=TRUE)
    genomepattern <- paste0(assembly, "$")
    bsgenome <- grep(genomepattern, available.genomes(), value=TRUE)
  }
  if (!bsgenome %in% installed.genomes()) {
    message("Your reads appear to be aligned against ", assembly)
    message("In order to correct for GC% and CpG density, you need a genome.")
    message("Our best guess for that genome is ", bsgenome)
    message("It doesn't look like you have installed it yet.")
    message("Try `BiocManager::install(\"", bsgenome, "\").)")
    stop("Need ", bsgenome, " to predict bin-level concentrations. Exiting.")
  }

  # Get read_count, fraglen, GC, and CpG_3 from GRanges object
  if (is.null(ssb_res)) ssb_res <- attr(fit, "data") # if feasible... ?
  reads <- subset(ssb_res$genomic, coverage != 0) # is this sensible (bias)?
  names(mcols(reads)) <- sub("coverage", "read_count", names(mcols(reads)))
  reads$fraglen <- width(reads)

  # can this be made superfluous?
  message("Adjusting for bin-level biases...")
  myGenome <- .load_genome(bsgenome)
  myBinSeqs <- Views(myGenome, reads)
  myBinGCs <- alphabetFrequency(myBinSeqs)[, c("G","C")]
  reads$GC <- rowSums(myBinGCs) / width(reads)
  myBinCpGs <- dinucleotideFrequency(myBinSeqs)[, "CG"]
  reads$CpG_3 <- myBinCpGs ** (1/3)  
  names(reads) <- as.character(reads) # coords 
  message("Done.")

  message("Starting pmol prediction...")
  reads$pred_conc <- predict(fit, as.data.frame(mcols(reads)))
  message("Done.")
 
  if (ret == "gr") return(reads)
  else return(.gr_prediction_to_df(reads, myBinSeqs))

} 


# helper fn
.load_genome <- function(bsgenome) {
  if (!any(grepl(bsgenome, loadedNamespaces()))) { 
    message("Attempting to load ", bsgenome, "...")
    x <- try(suppressMessages(attachNamespace(bsgenome)), silent=TRUE)
    if (inherits(x, "try-error")) stop("Unable to load library ", bsgenome, "!")
    message("OK.")
  }
  g <- try(getBSgenome(bsgenome, masked=FALSE, load.only=FALSE))
  if (inherits(g, "try-error")) stop("Unable to load ", bsgenome, ", exiting")
  return(g)
}


# helper fn
.gr_prediction_to_df <- function(gr, seqs) {

  res <- data.frame(chrom=decode(seqnames(gr)),
                    range.start=start(gr),
                    range.end=end(gr),
                    sequence=DNAStringSet(seqs),
                    read_count=gr$read_count,
                    fraglen=gr$fraglen,
                    GC=gr$GC,
                    CpG_3=gr$CpG_3,
                    pred_conc=gr$pred_conc)
  return(res)

}
