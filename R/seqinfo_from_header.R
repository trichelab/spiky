#' create seqinfo (and thus a standard chromosome filter) from a BAM/CRAM header
#' 
#' @param   x     the BAM/CRAM file or its header 
#' @param   gen   genome of the BAM/CRAM file, if known (NULL; autodetect)
#' @param   std   standard chromosomes only? (FALSE; will be empty if spikes)
#' @param   ret   return Seqinfo ("si", the default) or GRanges ("gr")? ("si")
#' 
#' @return        Seqinfo object or GRanges (or `as(seqinfo, "GRanges")`)
#'
#' @details 
#' 
#'   Setting std=TRUE on a spike-in CRAM or BAM will produce an empty result. 
#' 
#' @examples 
#' 
#'   library(Rsamtools)
#'   fl <- system.file("extdata", "ex1.bam", package="Rsamtools", mustWork=TRUE)
#'
#'   hdr <- scanBamHeader(BamFile(fl))
#'   si <- seqinfo_from_header(hdr)
#'   gr <- seqinfo_from_header(fl, ret="gr")
#'   stopifnot(identical(gr, as(si, "GRanges")))
#' 
#'   std_si <- seqinfo_from_header(fl, std=TRUE)
#'   seqlevels(std_si)
#'
#'   # for comparison with below
#'   data(spike, package="spiky") 
#'   spike 
#'
#'   # FIXME: add an example with a CRAM instead 
#'   sp <- system.file("extdata", "example.spike.bam", package="spiky")
#'   sp_gr <- seqinfo_from_header(sp, ret="gr")
#'   sp_gr 
#' 
#' @import Rsamtools
#' @import GenomeInfoDb 
#' 
#' @export 
seqinfo_from_header <- function(x, gen=NA, std=FALSE, ret=c("si", "gr")) {
  
  ret <- match.arg(ret)
  
  if (is.list(x) & "targets" %in% names(x)) { # it's a header, just convert 
    res <- .seqinfo_from_targets(x$targets, genome=gen) 
  } else {
    if (!is(x, "BamFile")) x <- BamFile(x)
    res <- .seqinfo_from_targets(scanBamHeader(x)$targets, genome=gen)
  }

  if (std) res <- keepStandardChromosomes(res)
  if (ret == "gr") res <- as(res, "GRanges")
  return(res)

}


# helper fn -- may break this out eventually
.seqinfo_from_targets <- function(chrs, genome=NA) { 

  if (is.null(names(chrs))) stop("Seqlengths have no contig names. Exiting.") 
  if (is.list(chrs) & "targets" %in% names(chrs)) chrs <- chrs$targets 
  res <- as(data.frame(chrom=names(chrs), seqlengths=chrs), "Seqinfo") 
 
  if (!is.na(genome) & !is.null(genome)) {

    genome(res) <- genome 

  } else { 

    # we can deduce a bit here, at least for mouse and human:
    if ("chrM" %in% seqlevels(res) | "MT" %in% seqlevels(res)) { 
      orig <- seqlevelsStyle(res)
      seqlevelsStyle(res) <- "UCSC"
      if ("chrM" %in% seqlevels(res)) {
        if (seqlengths(res)["chrM"] == 16299) genome(res) <- "mm10" 
        if (seqlengths(res)["chrM"] == 16571) genome(res) <- "hg19" # boooo
        if (seqlengths(res)["chrM"] == 16569) { 
          if ("chr1" %in% seqlevels(res) & 
              seqlengths(res)["chr1"] == 248956422) {
            genome(res) <- "hg38"
          } else { 
            genome(res) <- "GRCh37"
          }
        }
      }
      seqlevelsStyle(res) <- orig
    }

  }

  return(res) 

}
