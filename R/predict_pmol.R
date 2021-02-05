#' predict picomoles of DNA from a fit and readcounts 
#' 
#' this seems to just wrap output from bin_pmol and model_glm_pmol... ? 
#' 
#' @param   fit result of model_glm_pmol 
#' @param   ssb_res   the data / new data 
#' 
#' @return  data, a dataframe with sequence, read count, fragment length, GC% CpG cube rooted, and concentration prediction.
#' 
#' @import tools
#' @import reshape2
#' @import stringr
#' @import BSgenome 
#' @import S4Vectors
#' 
#' @examples 
#' 
#' data(ssb_res)
#' fit <- model_glm_pmol(covg_to_df(ssb_res)) 
#' pred <- predict_pmol(fit,ssb_res) 
#' 
#' @export 
predict_pmol <- function(fit, ssb_res=NULL, genome="BSgenome.Hsapiens.UCSC.hg38") {

  myGenome <- load_genome(genome)
  
  # Get read_count, fraglen, GC, and CpG_3 from GRanges object
  reads <-subset(ssb_res$genomic,coverage != 0)
  seqs <- getSeq(myGenome,reads)
  data <- data.frame("chrom" = decode(seqnames(reads)),"range" = decode(ranges(reads)),
                     "sequence" = seqs,"read_count" = reads$coverage)
  data$range.names = NULL
  data$range.width = NULL
  data$fraglen = str_length(data$sequence)
  data$GC = sapply(data$sequence,str_count,pattern="[CGcg]") / str_length(data$sequence)
  data$CpG_3 = sapply(data$sequence,str_count,pattern="CG") ^ (1/3) 
  
  if (is.null(ssb_res)) ssb_res <- attr(fit, "data") 
  message("start pmol prediction")
  data$pred_conc <- predict(fit, data)
  return(data) 

} 

load_genome <- function(genome) {
  try(suppressMessages(attachNamespace(genome)), silent=TRUE)
  myGenome <- getBSgenome(genome, masked=FALSE, load.only=FALSE)
  return(myGenome)
}
